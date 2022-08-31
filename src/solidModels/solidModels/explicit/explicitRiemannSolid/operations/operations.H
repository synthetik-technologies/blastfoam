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

Class
    Foam::operations

Description
    Define commonly used mathematical operations.

SourceFiles
    operations.C

\*---------------------------------------------------------------------------*/

#ifndef operations_H
#define operations_H

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class operations
{
    // Private data

        //- Mesh
         const fvMesh& mesh_;

        //- Eigen-values
        mutable vector eigVal_;

        //- Eigen-vectors
        mutable tensor eigVec_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        operations(const operations&);

        //- Disallow default bitwise assignment
        void operator=(const operations&);

public:

    //- Runtime type information
    TypeName("operations");

    // Constructors

        //- Construct from mesh
        operations
        (
            const fvMesh&
        );


    //- Destructor
    virtual ~operations();


    // Member functions

        //Access

            const vector& eigenValue() const
            {
                return eigVal_;
            }

            const tensor& eigenVector() const
            {
                return eigVec_;
            }

        // Edit

            //- Tensor cross product of two tensors
            tensor tensorProduct(const tensor&, const tensor&) const;

            // Compute eigen-values and eigen-vectors
            void eigenStructure(const tensor&) const;

            //- Inverse transpose of a volTensorField
            template<template<class> class Patch, class Mesh>
            tmp<GeometricField<tensor, Patch, Mesh>> invT
            (
                const GeometricField<tensor, Patch, Mesh>&
            ) const;

            //- Tensor cross product of two volTensorFields
            template<template<class> class Patch, class Mesh>
            tmp<GeometricField<tensor, Patch, Mesh>> tensorProduct
            (
                const GeometricField<tensor, Patch, Mesh>&,
                const GeometricField<tensor, Patch, Mesh>&
            ) const;

            //- Decompose a tensor field
            void decomposeTensor
            (
                const Field<tensor>&,
                Field<vector>&,
                Field<vector>&,
                Field<vector>&
            ) const;

            //- Decompose a geometric tensor field
            template<template<class> class Patch, class Mesh>
            void decomposeTensor
            (
                const GeometricField<tensor, Patch, Mesh>&,
                GeometricField<vector, Patch, Mesh>&,
                GeometricField<vector, Patch, Mesh>&,
                GeometricField<vector, Patch, Mesh>&
            ) const;

            //- Decompose tensor and return x row
            template<template<class> class Patch, class Mesh>
            tmp<GeometricField<vector, Patch, Mesh>> decomposeTensorX
            (
                const GeometricField<tensor, Patch, Mesh>&
            ) const;

            //- Decompose tensor and return y row
            template<template<class> class Patch, class Mesh>
            tmp<GeometricField<vector, Patch, Mesh>> decomposeTensorY
            (
                const GeometricField<tensor, Patch, Mesh>&
            ) const;

            //- Decompose tensor and return z row
            template<template<class> class Patch, class Mesh>
            tmp<GeometricField<vector, Patch, Mesh>> decomposeTensorZ
            (
                const GeometricField<tensor, Patch, Mesh>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Stabilised tensor inverse
tensor stabInv(const tensor&);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "operationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
