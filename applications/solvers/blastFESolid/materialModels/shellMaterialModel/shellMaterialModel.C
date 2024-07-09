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

#include "shellMaterialModel.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shellMaterialModel, 0);
    defineRunTimeSelectionTable(shellMaterialModel, linear);
    defineRunTimeSelectionTable(shellMaterialModel, nonLinear);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shellMaterialModel::shellMaterialModel
(
    const dictionary& dict,
    const feMesh1& mesh,
    const Field<vector>& D,
    const Field<vector>& U,
    const Field<vector>& theta,
    const Field<vector>& omega,
    const GeoType geoType
)
:
    materialModel(dict, mesh, D, U, false, geoType),
    theta_(theta),
    omega_(omega),
    h_(dict.lookup<scalar>("thickness"))
{
    // Reset zone id to use face zones
    id_ = mesh.mesh().faceZones().findZoneID(name_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shellMaterialModel::~shellMaterialModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::vector, 3> Foam::shellMaterialModel::covariantBasis
(
    const tensor& gradX,
    const vector& v
) const
{
    return FixedList<vector, 3>({gradX.x(), gradX.y(), v});
}


Foam::FixedList<Foam::vector, 3>
Foam::shellMaterialModel::contravariantBasis
(
    const FixedList<vector, 3>& cov
) const
{
    const scalar denom(1.0/(cov[0] & (cov[1] ^ cov[2])));

    return FixedList<vector, 3>
    (
        {
            (cov[1] ^ cov[2])*denom,
            (cov[2] ^ cov[0])*denom,
            (cov[0] ^ cov[1])*denom
        }
    );
}


Foam::FixedList<Foam::vector, 3>
Foam::shellMaterialModel::localBasis
(
    const FixedList<vector, 3>& cov
) const
{
    const vector e3 = normalised(cov[2]);
    const vector e1 = normalised(cov[1] ^ e3);
    return FixedList<vector, 3>({e1, e3 ^ e1, e3});
}


const Foam::labelList& Foam::shellMaterialModel::elements() const
{
    return
        id_ >= 0
      ? static_cast<const labelList&>(mesh_.mesh().faceZones()[id_])
      : labelList::null();
}

void Foam::shellMaterialModel::updateStrain
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


Foam::scalar Foam::shellMaterialModel::addForce
(
    UList<vector>& force,
    UList<vector>& moment,
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
    UIndirectList<vector> moment_loc(moment, elem);
    UIndirectList<symmTensor> sigma_loc(sigma, elem);

    scalar vp = 0.0;
    symmTensor s;
    if (geoType_ == TOTAL_LAGRANGIAN)
    {
        tensor P;
        forAll(elem.ir(), rulei)
        {
            const vector n(elem.calcOrtho(mesh_.Js()[elemi][rulei]));
            const integrationPoint& ip = elem.ir()[rulei];

            // Compute stress at the integration point
            scalar vp_i = this->calcPiola(P, s, elemi, rulei);

            const scalar w = ip.w()*Ws[rulei];
            const scalarList& shape = shapes[rulei];
            const scalarRectangularMatrix& B = Bs[rulei];
            const scalar tH = h_*elem.ir()[rulei].z()/2.0;
            for (label si = 0; si < B.m(); si++)
            {
                vp += shape[si]*vp_i;
                sigma_loc[si] += w*shape[si]*s;

                force_loc[si] -=
                    w*mag(tH)
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

                // Shear stress
                const tensor tau(P - (P & n)*n);
                moment_loc[si] -=
                    w*sqr(tH)
                   *vector
                    (
                        B(si, 0)*tau.xx()
                      + B(si, 1)*tau.xy()
                      + B(si, 2)*tau.xz(),

                        B(si, 0)*tau.yx()
                      + B(si, 1)*tau.yy()
                      + B(si, 2)*tau.yz(),

                        B(si, 0)*tau.zx()
                      + B(si, 1)*tau.zy()
                      + B(si, 2)*tau.zz()
                    );
            }
        }
    }
    else
    {
        forAll(elem.ir(), rulei)
        {
            const vector n(elem.calcOrtho(mesh_.Js()[elemi][rulei]));
            const integrationPoint& ip = elem.ir()[rulei];

            // Compute stress at the integration point
            scalar vp_i = this->calcSigma(s, elemi, rulei);

            const scalar w = ip.w()*Ws[rulei];
            const scalarList& shape = shapes[rulei];
            const scalarRectangularMatrix& B = Bs[rulei];
            const scalar tH = h_*elem.ir()[rulei].z()/2.0;
            for (label si = 0; si < B.m(); si++)
            {
                vp += shape[si]*vp_i;
                sigma_loc[si] += w*shape[si]*s;
                force_loc[si] -=
                    w*mag(tH)
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

                // Shear stress
                const tensor tau(s - (s & n)*n);
                moment_loc[si] -=
                    w*sqr(tH)
                   *vector
                    (
                        B(si, 0)*tau.xx()
                      + B(si, 1)*tau.xy()
                      + B(si, 2)*tau.xz(),

                        B(si, 0)*tau.xy()
                      + B(si, 1)*tau.yy()
                      + B(si, 2)*tau.yz(),

                        B(si, 0)*tau.zx()
                      + B(si, 1)*tau.zy()
                      + B(si, 2)*tau.zz()
                    );
            }
        }
    }

    update_ = false;
    return vp;
}


// ************************************************************************* //

