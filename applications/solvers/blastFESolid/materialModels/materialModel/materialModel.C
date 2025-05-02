/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "materialModel.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materialModel, 0);
    defineRunTimeSelectionTable(materialModel, linear);
    defineRunTimeSelectionTable(materialModel, nonLinear);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::materialModel::Kmu_Enu
(
    const scalar E, const scalar nu,
    scalar& K, scalar& mu,
    const bool ps
)
{
    mu = E/(2.0*(1.0 + nu));
    if (ps)
    {
        K = E/(2.0*(1.0 - nu));
//         K = nu*E/((1.0 + nu)*(1.0 - nu)) + 2.0*mu/3.0;
    }
    else
    {
        K = E/(3.0*(1.0 - 2.0*nu));
//         K = nu*E/((1.0 + nu)*(1.0 - 2.0*nu)) + 2.0*mu/3.0;
    }
}


void Foam::materialModel::Enu_Kmu
(
    const scalar K, const scalar mu,
    scalar& E, scalar& nu,
    const bool ps
)
{
    if (ps)
    {
        E = 4.0*K*mu/(K + mu);
        nu = (K - mu)/(K + mu);
    }
    else
    {
        E = 9.0*K*mu/(3.0*K + mu);
        nu = (3.0*K - 2.0*mu)/(6.0*K + 2.0*mu);
    }
}

Foam::scalar Foam::materialModel::lambda_Kmu
(
    const scalar K, const scalar mu,
    const bool ps
)
{
    if (ps)
    {
        return K - mu;
    }
    else
    {
        return K - 2.0*mu/3.0;
    }
}

Foam::scalar Foam::materialModel::lambda_Enu
(
    const scalar E, const scalar nu,
    const bool ps
)
{
    if (ps)
    {
        return E*nu/((1.0 + nu)*(1.0 - nu));
    }
    else
    {
        return E*nu/((1.0 + nu)*(1.0 - 2.0*nu));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModel::materialModel
(
    const dictionary& dict,
    const feMesh1& mesh,
    const pointVectorField& D,
    const pointVectorField& U,
    const bool planeStress,
    const GeoType geoType
)
:
    mesh_(mesh),
    geoType_(geoType),
    D_(D),
    U_(U),
    rho_("rho", dimDensity, 0.0),
    name_(dict.name()),
    planeStress_(planeStress),
    update_(false)
{
    rho_.read(dict);

    if (planeStress_)
    {
        const polyMesh& pmesh = mesh.mesh();
        const Vector<label>& solutionD = pmesh.solutionD();
        label nD = 0;
        if (solutionD[0] < 0)
        {
            nD++;
            planeStressDir_ = symmTensor::XX;
            nonPlaneStressDir1_ = direction(symmTensor::YY);
            nonPlaneStressDir2_ = direction(symmTensor::ZZ);
        }
        if (solutionD[1] < 0)
        {
            nD++;
            nonPlaneStressDir1_ = direction(symmTensor::XX);
            planeStressDir_ = symmTensor::YY;
            nonPlaneStressDir2_ = direction(symmTensor::ZZ);
        }
        if (solutionD[2] < 0)
        {
            nD++;
            nonPlaneStressDir1_ = direction(symmTensor::XX);
            nonPlaneStressDir2_ = direction(symmTensor::YY);
            planeStressDir_ = symmTensor::ZZ;
        }
        if (nD != 1)
        {
            FatalErrorInFunction
                << "For planeStress, this material law assumes one empty "
                << "direction, but " << nD << "  empty directions were found."
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModel::~materialModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::materialModel::elements() const
{
    return
        mesh_.mesh().cellZones().found(name_)
      ? static_cast<const labelList&>(mesh_.mesh().cellZones()[name_])
      : labelList::null();
}


void Foam::materialModel::updateStrain
(
    symmTensor& epsilon,
    const scalar coeff,
    const symmTensor& sigma
) const
{
    if (planeStress_)
    {
        epsilon[planeStressDir_] =
          - coeff
           *(
                sigma[nonPlaneStressDir1_]
              + sigma[nonPlaneStressDir2_]
            );
    }
}


Foam::scalar Foam::materialModel::addForce
(
    UList<vector>& force,
    UList<symmTensor>& sigma,
    const label elemi
)
{
    update_ = true;

    const element& elem = mesh_.elements()[elemi];
    const List<scalarList>& shapes = mesh_.shapes()[elemi];
    const List<scalar>& Ws = mesh_.Ws()[elemi];
    const List<scalarRectangularMatrix>& Bs = mesh_.Bs()[elemi];
    UIndirectList<vector> force_loc(force, elem);
    UIndirectList<symmTensor> sigma_loc(sigma, elem);

    scalar vp = 0.0;
    symmTensor s;
    if (geoType_ == TOTAL_LAGRANGIAN)
    {
        tensor P;
        forAll(elem.ir(), rulei)
        {
            const integrationPoint& ip = elem.ir()[rulei];

            // Compute stress at the integration point
            scalar vp_i = this->calcPiola(P, s, elemi, rulei);

            const scalar w = ip.w()*Ws[rulei];
            const scalarList& shape = shapes[rulei];
            const scalarRectangularMatrix& B = Bs[rulei];
            for (label si = 0; si < B.m(); si++)
            {
                vp += shape[si]*vp_i;
                sigma_loc[si] += w*shape[si]*s;
                force_loc[si] -=
                    w
                   *vector
                    (
                        B(si, 0)*P.xx()
                      + B(si, 1)*P.xy()
                      + B(si, 2)*P.xz(),

                        B(si, 0)*P.yx()
                      + B(si, 1)*P.yy()
                      + B(si, 2)*P.yz(),

                        B(si, 0)*P.zx()
                      + B(si, 1)*P.zy()
                      + B(si, 2)*P.zz()
                    );
            }
        }
    }
    else
    {
        forAll(elem.ir(), rulei)
        {
            const integrationPoint& ip = elem.ir()[rulei];

            // Compute stress at the integration point
            scalar vp_i = this->calcSigma(s, elemi, rulei);

            const scalar w = ip.w()*Ws[rulei];
            const scalarList& shape = shapes[rulei];
            const scalarRectangularMatrix& B = Bs[rulei];
            for (label si = 0; si < B.m(); si++)
            {
                vp += shape[si]*vp_i;
                sigma_loc[si] += w*shape[si]*s;
                force_loc[si] -=
                    w
                   *vector
                    (
                        B(si, 0)*s.xx()
                      + B(si, 1)*s.xy()
                      + B(si, 2)*s.xz(),

                        B(si, 0)*s.xy()
                      + B(si, 1)*s.yy()
                      + B(si, 2)*s.yz(),

                        B(si, 0)*s.xz()
                      + B(si, 1)*s.yz()
                      + B(si, 2)*s.zz()
                    );
            }
        }
    }

    update_ = false;
    return vp;
}


// ************************************************************************* //

