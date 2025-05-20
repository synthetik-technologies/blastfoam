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

    viscousPressure_(dict.lookup<bool>("viscousPressure")),
    linearBulkViscosityCoeff_
    (
        dimensionedScalar::lookupOrDefault
        (
            "linearBulkViscosityCoeff",
            dict,
            0.06
        )
    ),
    quadraticBulkViscosityCoeff_
    (
        dimensionedScalar::lookupOrDefault
        (
            "quadraticBulkViscosityCoeff",
            dict,
            1.2
        )
    ),

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
    const label elemi,
    const scalar L
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

            if (viscousPressure_)
            {
                const scalar visP = this->viscousPressure
                (
                    vp_i,
                    L,
                    elemi,
                    rulei
                );
                P.xx() += visP;
                P.yy() += visP;
                P.zz() += visP;
            }

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

            if (viscousPressure_)
            {
                const scalar visP = this->viscousPressure
                (
                    vp_i,
                    L,
                    elemi,
                    rulei
                );
                s.xx() += visP;
                s.yy() += visP;
                s.zz() += visP;
            }

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

Foam::scalar Foam::materialModel::viscousPressure
(
    const scalar vp,
    const scalar L,
    const label elemi,
    const label rulei
) const
{
    // Element displacment and velocity
    const UIndirectList<vector> U(this->U_, mesh_.elements()[elemi]);

    // Element data
    const scalarRectangularMatrix& dshape = mesh_.dshapes()[elemi][rulei];
    const tensor& invJ = mesh_.invJs()[elemi][rulei];
    const scalar epsilonDot(tr(gradient(U, dshape, invJ))/3.0);
    return
        rho_.value()
       *(
            linearBulkViscosityCoeff_.value()*epsilonDot*vp*L
          + sqr(max(epsilonDot, 0.0)*quadraticBulkViscosityCoeff_.value()*L)
        );

}


void Foam::materialModel::eigenStructure
(
    const tensor& ten,
    vector& eigVal,
    tensor& eigVec
)
{
    tensor t(ten);

    scalar it_max = 100;
    scalar it_num = 0;
    scalar rot_num = 0;
    scalar size = 3;
    scalar gapq = 0.0;

    label i,j;
    label k,l,m,p1,q = 0;

    scalar term, termp, termq = 0.0;
    scalar g,h,c,w1 = 0.0;
    scalar theta,thresh = 0.0;
    scalar s,t1,tau1 = 0.0;

    vector d = vector(t.xx(),t.yy(),t.zz());
    vector bw = vector(t.xx(),t.yy(),t.zz());
    vector zw = vector::zero;
    tensor v1 = tensor::I;

    while (it_num++ < it_max)
    {
        // The convergence threshold is based on the size of the elements in
        // the strict upper triangle of the matrix.
        thresh = 0.0;

        for (j = 0; j < size; j++)
        {
            for (i = 0; i < j; i++)
            {
                thresh += sqr(t[i + j*size]);
            }
        }

        thresh = sqrt(thresh)/(4.0*size);

        if ( thresh < small )
        {
          break;
        }

        for (p1 = 0; p1 < size; p1++)
        {
            for (q = p1 + 1; q < size; q++)
            {
                gapq = 10.0*mag(t[p1 + q*size]);
                termp = gapq + mag(d[p1]);
                termq = gapq + mag(d[q]);

                // Annihilate tiny off-diagonal elements
                if (4 < it_num && termp == mag(d[p1]) && termq == mag(d[q]))
                {
                  t[p1+q*size] = 0.0;
                }

                //  Otherwise, apply a rotation
                else if (thresh <= mag(t[p1 + q*size]))
                {
                    h = d[q] - d[p1];
                    term = mag(h) + gapq;

                    if (term == mag(h))
                    {
                        t1 = t[p1 + q*size]/h;
                    }
                    else
                    {
                        theta = 0.5 * h/t[p1 + q*size];
                        t1 = 1.0/(mag(theta) + sqrt(1.0 + theta*theta));

                        if (theta < 0.0)
                        {
                            t1 = - t1;
                        }
                    }

                    c = 1.0/sqrt(1.0 + t1*t1);
                    s = t1*c;
                    tau1 = s/(1.0 + c);
                    h = t1*t[p1 + q*size];

                    //  Accumulate corrections to diagonal elements
                    zw[p1] = zw[p1] - h;
                    zw[q] = zw[q] + h;
                    d[p1] = d[p1] - h;
                    d[q] = d[q] + h;
                    t[p1 + q*size] = 0.0;

                    // Rotate, using information from the upper triangle of A
                    // only
                    for (j = 0; j < p1; j++)
                    {
                        g = t[j + p1*size];
                        h = t[j + q*size];
                        t[j + p1*size] = g - s*(h + g*tau1);
                        t[j + q*size] = h + s*(g - h*tau1);
                    }

                    for (j = p1 + 1; j < q; j++)
                    {
                        g = t[p1 + j*size];
                        h = t[j + q*size];
                        t[p1 + j*size] = g - s*(h + g*tau1);
                        t[j + q*size] = h + s*(g - h*tau1);
                    }

                    for (j = q + 1; j < size; j++)
                    {
                        g = t[p1 + j*size];
                        h = t[q + j*size];
                        t[p1 + j*size] = g - s*(h + g*tau1);
                        t[q + j*size] = h + s*(g - h*tau1);
                    }

                    //  Accumulate information in the eigenvector matrix
                    for (j = 0; j < size; j++)
                    {
                        g = v1[j + p1*size];
                        h = v1[j + q*size];
                        v1[j + p1*size] = g - s*(h + g*tau1);
                        v1[j + q*size] = h + s*(g - h*tau1);
                    }
                    rot_num = rot_num + 1;
                }
            }
        }
    }

    // Restore upper triangle of input matrix
    for (i = 0; i < size; i++)
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }

    //  Ascending sort the eigenvalues and eigenvectors
    for (k = 0; k < size - 1; k++)
    {
        m = k;
        for (l = k + 1; l < size; l++)
        {
            if (d[l] < d[m])
            {
                m = l;
            }
        }

        if (m != k)
        {
            t1 = d[m];
            d[m] = d[k];
            d[k] = t1;
            for (i = 0; i < size; i++)
            {
                w1 = v1[i + m*size];
                v1[i + m*size] = v1[i + k*size];
                v1[i + k*size] = w1;
            }
        }
    }

    // Corrections for calculating inverse
    tensor sub(tensor::zero);

    for (i = 0; i < 3; i++)
    {
        if (d[i] < SMALL)
        {
            d[i] = 1;
            sub +=
                d[i]
               *vector(v1[3*i], v1[3*i + 1], v1[3*i + 2])
               *vector(v1[3*i], v1[3*i + 1], v1[3*i + 2]);
        }
    }

    eigVal = d;
    eigVec = v1;
}



// ************************************************************************* //

