/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam::fvc

Description
    Calculate the surface gradient of the given field.

SourceFiles
    fvcGradf.C

\*---------------------------------------------------------------------------*/


#ifndef fvcGradf_H
#define fvcGradf_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvsPatchField,
            surfaceMesh
        >
    > fGrad
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const GeometricField<Type, pointPatchField, pointMesh>& pf,
        bool interpolate = false
    );

//     template<class Type>
//     tmp
//     <
//         GeometricField
//         <
//             typename outerProduct<vector, Type>::type,
//             fvsPatchField,
//             surfaceMesh
//         >
//     > fGrad
//     (
//         const GeometricField<Type, fvPatchField, volMesh>& vf,
//         const GeometricField<Type, pointPatchField, pointMesh>& pf,
//         const GeometricField
//         <
//             typename outerProduct<vector, Type>::type,
//             fvPatchField,
//             volMesh
//         >& gradVf
//     );

    template<class Type>
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvsPatchField,
            surfaceMesh
        >
    > fsGrad
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const GeometricField<Type, pointPatchField, pointMesh>& pf
    );

    template<class Type, template<class> class FaceList>
    tmp<Field<typename outerProduct<vector, Type>::type> > fGrad
    (
        const PrimitivePatch<FaceList<face>, const pointField&>& patch,
        const Field<Type>& ppf
    );

    template<class Type>
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    > grad
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const GeometricField<Type, pointPatchField, pointMesh>& pf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvcGradf.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
