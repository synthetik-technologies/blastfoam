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


tmp<volTensorField> operations::invT
(
    const GeometricField<tensor, fvPatchField, volMesh>& t
) const
{
    return Foam::T(Foam::inv(t));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> operations::tensorProduct
(
    const GeometricField<tensor, fvPatchField, volMesh>& T1,
    const GeometricField<tensor, fvPatchField, volMesh>& T2
) const
{
    tmp<GeometricField<tensor, fvPatchField, volMesh> > tsf
    (
        GeometricField<tensor, fvPatchField, volMesh>::New
        (
            "P",
            mesh_,
            dimensioned<tensor>
            (
                "P",
                T1.dimensions()*T2.dimensions(),
                pTraits<tensor>::one
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh>& P = tsf.ref();

    forAll(mesh_.cells(), cellID)
    {
        P[cellID] = operations::tensorProduct(T1[cellID], T2[cellID]);
    }

    return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    const GeometricField<tensor, fvPatchField, volMesh>& T,
    GeometricField<vector, fvPatchField, volMesh>& Vx,
    GeometricField<vector, fvPatchField, volMesh>& Vy,
    GeometricField<vector, fvPatchField, volMesh>& Vz
) const
{
    forAll(mesh_.cells(), cellID)
    {
        Vx[cellID] = vector(T[cellID].xx(), T[cellID].xy(), T[cellID].xz());
        Vy[cellID] = vector(T[cellID].yx(), T[cellID].yy(), T[cellID].yz());
        Vz[cellID] = vector(T[cellID].zx(), T[cellID].zy(), T[cellID].zz());
    }

    if (Pstream::parRun())
    {
        Vx.correctBoundaryConditions();
        Vy.correctBoundaryConditions();
        Vz.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void operations::decomposeTensor
(
    const GeometricField<tensor, fvsPatchField, surfaceMesh>& T,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Vx,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Vy,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Vz
) const
{
    forAll(mesh_.owner(), faceID)
    {
        Vx[faceID] = vector(T[faceID].xx(), T[faceID].xy(), T[faceID].xz());
        Vy[faceID] = vector(T[faceID].yx(), T[faceID].yy(), T[faceID].yz());
        Vz[faceID] = vector(T[faceID].zx(), T[faceID].zy(), T[faceID].zz());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> operations::decomposeTensorX
(
    const GeometricField<tensor, fvPatchField, volMesh>& T
) const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            T.name() + ".x",
            mesh_,
            dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& Vx = tvf.ref();

    forAll(mesh_.cells(), cellID)
    {
        Vx[cellID] = vector(T[cellID].xx(), T[cellID].xy(), T[cellID].xz());
    }

    if (Pstream::parRun())
    {
        Vx.correctBoundaryConditions();
    }

    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> operations::decomposeTensorY
(
    const GeometricField<tensor, fvPatchField, volMesh>& T
) const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            T.name() + ".y",
            mesh_,
            dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& Vy = tvf.ref();

    forAll(mesh_.cells(), cellID)
    {
        Vy[cellID] = vector(T[cellID].yx(), T[cellID].yy(), T[cellID].yz());
    }

    if (Pstream::parRun())
    {
        Vy.correctBoundaryConditions();
    }

    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> operations::decomposeTensorZ
(
    const GeometricField<tensor, fvPatchField, volMesh>& T
) const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            T.name() + ".z",
            mesh_,
            dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& Vz = tvf.ref();

    forAll(mesh_.cells(), cellID)
    {
        Vz[cellID] = vector(T[cellID].zx(), T[cellID].zy(), T[cellID].zz());
    }

    if (Pstream::parRun())
    {
        Vz.correctBoundaryConditions();
    }

    return tvf;
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

    while (it_num < it_max)
    {
        it_num += 1;

        // The convergence threshold is based on the size of the elements in
        // the strict upper triangle of the matrix.
        thresh = 0.0;

        for (j = 0; j < size; j++)
        {
            for (i = 0; i < j; i++)
            {
                thresh = thresh + t[i + j*size] * t[i + j*size];
            }
        }

        thresh = ::sqrt(thresh)/(4*size);

        if ( thresh == 0.0 )
        {
          break;
        }

        for (p1 = 0; p1 < size; p1++)
        {
            for (q = p1 + 1; q < size; q++)
            {
                gapq = 10.0 * fabs(t[p1 + q*size]);
                termp = gapq + fabs(d[p1]);
                termq = gapq + fabs(d[q]);

                // Annihilate tiny off-diagonal elements
                if (4 < it_num && termp == fabs(d[p1]) && termq == fabs(d[q]))
                {
                  t[p1+q*size] = 0.0;
                }

                //  Otherwise, apply a rotation
                else if (thresh <= fabs(t[p1 + q*size]))
                {
                    h = d[q] - d[p1];
                    term = fabs(h) + gapq;

                    if (term == fabs(h))
                    {
                        t1 = t[p1 + q*size]/h;
                    }
                    else
                    {
                        theta = 0.5 * h/t[p1 + q*size];
                        t1 = 1.0/(fabs(theta) + ::sqrt(1.0 + theta*theta));

                        if (theta < 0.0)
                        {
                            t1 = - t1;
                        }
                    }

                    c = 1.0/::sqrt(1.0 + t1*t1);
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

    for (i=0; i<3; i++)
    {
        if (d[i] < SMALL)
        {
            d[i] = 1;
            sub +=
                d[i] * vector(v1[3*i],v1[3*i+1],v1[3*i+2])
               *vector(v1[3*i],v1[3*i+1],v1[3*i+2]);
        }
    }

    eigVal_ = d;
    eigVec_ = v1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
