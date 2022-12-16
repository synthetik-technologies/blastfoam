/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "unsExplicitNonLinearSolid.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

void Foam::solidModels::unsExplicitNonLinearSolid::updateWavespeeds()
{
    wavespeed_ =
        fvc::interpolate(sqrt(this->mechanical().elasticModulus()/this->rho()));
    sWavespeed_ =
        fvc::interpolate
        (
            sqrt(this->mechanical().shearModulus()/this->rho())
        );

    vector eigVal;
    tensor eigVec;

    surfaceTensorField Ff(this->Ff());
    surfaceTensorField Cf(Ff.T() & Ff);
    forAll(Cf, facei)
    {
        eigenStructure(Cf[facei], eigVal, eigVec);
        scalar s = sqrt(cmptMin(eigVal));
        this->wavespeed_[facei] /= s;
        this->sWavespeed_[facei] /= s;
    }

    const surfaceTensorField::Boundary& bCf(Cf.boundaryField());
    surfaceScalarField::Boundary& bwavespeed =
        this->wavespeed_.boundaryFieldRef();
    surfaceScalarField::Boundary& bsWavespeed =
        this->sWavespeed_.boundaryFieldRef();
    forAll(bCf, patchi)
    {
        forAll(bCf[patchi], facei)
        {
            eigenStructure(bCf[patchi][facei], eigVal, eigVec);
            scalar s = sqrt(cmptMin(eigVal));
            bwavespeed[patchi][facei] /= s;
            bsWavespeed[patchi][facei] /= s;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModels::unsExplicitNonLinearSolid::unsExplicitNonLinearSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    unsExplicitSolid(type, mesh, nonLinear, isSolid)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::solidModels::unsExplicitNonLinearSolid::eigenStructure
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
