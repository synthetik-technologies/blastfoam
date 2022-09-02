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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<template<class> class Patch, class Mesh>
tmp<GeometricField<tensor, Patch, Mesh>> operations::invT
(
    const GeometricField<tensor, Patch, Mesh>& t
) const
{
    return Foam::T(Foam::inv(t));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Patch, class Mesh>
tmp<GeometricField<tensor, Patch, Mesh>> operations::tensorProduct
(
    const GeometricField<tensor, Patch, Mesh>& T1,
    const GeometricField<tensor, Patch, Mesh>& T2
) const
{
    tmp<GeometricField<tensor, Patch, Mesh>> tP
    (
        GeometricField<tensor, Patch, Mesh>::New
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
    GeometricField<tensor, Patch, Mesh>& P = tP.ref();

    forAll(T1, i)
    {
        P[i] = tensorProduct(T1[i], T2[i]);
    }

    typename GeometricField<tensor, Patch, Mesh>::Boundary& pP =
        P.boundaryFieldRef();
    forAll(pP, patchi)
    {
        forAll(pP[patchi], facei)
        {
            pP[patchi][facei] =
                tensorProduct
                (
                    T1.boundaryField()[patchi][facei],
                    T2.boundaryField()[patchi][facei]
                );
        }
    }

    return tP;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Patch, class Mesh>
void operations::decomposeTensor
(
    const GeometricField<tensor, Patch, Mesh>& T,
    GeometricField<vector, Patch, Mesh>& Tx,
    GeometricField<vector, Patch, Mesh>& Ty,
    GeometricField<vector, Patch, Mesh>& Tz
) const
{
    decomposeTensor
    (
        T.primitiveField(),
        Tx.primitiveFieldRef(),
        Ty.primitiveFieldRef(),
        Tz.primitiveFieldRef()
    );
    const typename GeometricField<tensor, Patch, Mesh>::Boundary& pT =
        T.boundaryField();
    typename GeometricField<vector, Patch, Mesh>::Boundary& pTx =
        Tx.boundaryFieldRef();
    typename GeometricField<vector, Patch, Mesh>::Boundary& pTy =
        Ty.boundaryFieldRef();
    typename GeometricField<vector, Patch, Mesh>::Boundary& pTz =
        Tz.boundaryFieldRef();
    forAll(pT, patchi)
    {
        decomposeTensor
        (
            pT[patchi],
            pTx[patchi],
            pTy[patchi],
            pTz[patchi]
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Patch, class Mesh>
tmp<GeometricField<vector, Patch, Mesh>> operations::decomposeTensorX
(
    const GeometricField<tensor, Patch, Mesh>& T
) const
{
    tmp<GeometricField<vector, Patch, Mesh>> tvf
    (
        GeometricField<vector, Patch, Mesh>::New
        (
            T.name() + ".x",
            mesh_,
            dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, Patch, Mesh>& Tx = tvf.ref();

    forAll(T, cellID)
    {
        Tx[cellID] = vector(T[cellID].xx(), T[cellID].xy(), T[cellID].xz());
    }

    const typename GeometricField<tensor, Patch, Mesh>::Boundary& pT =
        T.boundaryField();
    typename GeometricField<vector, Patch, Mesh>::Boundary& pTx =
        Tx.boundaryFieldRef();
    forAll(pT, patchi)
    {
        forAll(pT[patchi], facei)
        {
            pTx[patchi][facei] =
                vector
                (
                    pT[patchi][facei].xx(),
                    pT[patchi][facei].xy(),
                    pT[patchi][facei].xz()
                );
        }
    }
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Patch, class Mesh>
tmp<GeometricField<vector, Patch, Mesh>> operations::decomposeTensorY
(
    const GeometricField<tensor, Patch, Mesh>& T
) const
{
    tmp<GeometricField<vector, Patch, Mesh>> tvf
    (
        GeometricField<vector, Patch, Mesh>::New
        (
            T.name() + ".y",
            mesh_,
            dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, Patch, Mesh>& Ty = tvf.ref();

    forAll(T, cellID)
    {
        Ty[cellID] = vector(T[cellID].yx(), T[cellID].yy(), T[cellID].yz());
    }

    const typename GeometricField<tensor, Patch, Mesh>::Boundary& pT =
        T.boundaryField();
    typename GeometricField<vector, Patch, Mesh>::Boundary& pTy =
        Ty.boundaryFieldRef();
    forAll(pT, patchi)
    {
        forAll(pT[patchi], facei)
        {
            pTy[patchi][facei] =
                vector
                (
                    pT[patchi][facei].yx(),
                    pT[patchi][facei].yy(),
                    pT[patchi][facei].yz()
                );
        }
    }
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Patch, class Mesh>
tmp<GeometricField<vector, Patch, Mesh>> operations::decomposeTensorZ
(
    const GeometricField<tensor, Patch, Mesh>& T
) const
{
    tmp<GeometricField<vector, Patch, Mesh>> tvf
    (
        GeometricField<vector, Patch, Mesh>::New
        (
            T.name() + ".z",
            mesh_,
            dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, Patch, Mesh>& Tz = tvf.ref();

    forAll(T, cellID)
    {
        Tz[cellID] = vector(T[cellID].zx(), T[cellID].zy(), T[cellID].zz());
    }

    const typename GeometricField<tensor, Patch, Mesh>::Boundary& pT =
        T.boundaryField();
    typename GeometricField<vector, Patch, Mesh>::Boundary& pTz =
        Tz.boundaryFieldRef();
    forAll(pT, patchi)
    {
        forAll(pT[patchi], facei)
        {
            pTz[patchi][facei] =
                vector
                (
                    pT[patchi][facei].zx(),
                    pT[patchi][facei].zy(),
                    pT[patchi][facei].zz()
                );
        }
    }
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
