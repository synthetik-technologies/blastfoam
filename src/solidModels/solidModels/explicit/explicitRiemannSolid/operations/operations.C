/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(operations, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

operations::operations
(
    const fvMesh& vm
)
:
    mesh_(vm)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

operations::~operations()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tensor operations::tensorProduct
(
    const tensor& T1, const tensor& T2
) const
{
    tensor P = tensor::zero;

    P.xx() =
        (T1.yy()*T2.zz()) - (T1.yz()*T2.zy()) + (T1.zz()*T2.yy())
      - (T1.zy()*T2.yz());

    P.xy() =
        (T1.yz()*T2.zx()) - (T1.yx()*T2.zz()) + (T1.zx()*T2.yz())
      - (T1.zz()*T2.yx());

    P.xz() =
        (T1.yx()*T2.zy()) - (T1.yy()*T2.zx()) + (T1.zy()*T2.yx())
      - (T1.zx()*T2.yy());

    P.yx() =
        (T1.xz()*T2.zy()) - (T1.xy()*T2.zz()) + (T1.zy()*T2.xz())
      - (T1.zz()*T2.xy());

    P.yy() =
        (T1.zz()*T2.xx()) - (T1.zx()*T2.xz()) + (T1.xx()*T2.zz())
      - (T1.xz()*T2.zx());

    P.yz() =
        (T1.zx()*T2.xy()) - (T1.zy()*T2.xx()) + (T1.xy()*T2.zx())
      - (T1.xx()*T2.zy());

    P.zx() =
        (T1.xy()*T2.yz()) - (T1.xz()*T2.yy()) + (T1.yz()*T2.xy())
      - (T1.yy()*T2.xz());

    P.zy() =
        (T1.xz()*T2.yx()) - (T1.xx()*T2.yz()) + (T1.yx()*T2.xz())
      - (T1.yz()*T2.xx());

    P.zz() =
        (T1.xx()*T2.yy()) - (T1.xy()*T2.yx()) + (T1.yy()*T2.xx())
      - (T1.yx()*T2.xy());

    return P;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void operations::decomposeTensor
(
    const Field<tensor>& T,
    Field<vector>& Tx,
    Field<vector>& Ty,
    Field<vector>& Tz
) const
{
    forAll(T, i)
    {
        Tx[i] = vector(T[i].xx(), T[i].xy(), T[i].xz());
        Ty[i] = vector(T[i].yx(), T[i].yy(), T[i].yz());
        Tz[i] = vector(T[i].zx(), T[i].zy(), T[i].zz());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void operations::eigenStructure(const tensor& ten) const
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

    eigVal_ = d;
    eigVec_ = v1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tensor stabInv(const tensor& t)
{
    if (magSqr(t) < small)
    {
        return tensor::zero;
    }

    scalar scale = magSqr(t);
    Vector<bool> removeCmpts
    (
        magSqr(t.xx())/scale < small,
        magSqr(t.yy())/scale < small,
        magSqr(t.zz())/scale < small
    );
    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        tensor tPlus(t);

        if (removeCmpts.x())
        {
            tPlus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tPlus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tPlus += tensor(0,0,0,0,0,0,0,0,1);
        }

        tensor tInv = inv(tPlus);

        if (removeCmpts.x())
        {
            tInv -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tInv -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tInv -= tensor(0,0,0,0,0,0,0,0,1);
        }
        return tInv;
    }
    return inv(t);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
