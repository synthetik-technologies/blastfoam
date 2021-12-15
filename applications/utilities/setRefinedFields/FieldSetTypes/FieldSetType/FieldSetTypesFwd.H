/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef FieldSetTypesFwd_H
#define FieldSetTypesFwd_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldSetType.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef FieldSetType<scalar, fvPatchField, volMesh> scalarVolFieldSetType;
typedef FieldSetType<vector, fvPatchField, volMesh> vectorVolFieldSetType;
typedef FieldSetType<symmTensor, fvPatchField, volMesh>
    symmTensorVolFieldSetType;
typedef FieldSetType<sphericalTensor, fvPatchField, volMesh>
    sphericalTensorVolFieldSetType;
typedef FieldSetType<tensor, fvPatchField, volMesh> tensorVolFieldSetType;

typedef FieldSetType<scalar, fvsPatchField, surfaceMesh>
    scalarSurfaceFieldSetType;
typedef FieldSetType<vector, fvsPatchField, surfaceMesh>
    vectorSurfaceFieldSetType;
typedef FieldSetType<symmTensor, fvsPatchField, surfaceMesh>
    symmTensorSurfaceFieldSetType;
typedef FieldSetType<sphericalTensor, fvsPatchField, surfaceMesh>
    sphericalTensorSurfaceFieldSetType;
typedef FieldSetType<tensor, fvsPatchField, surfaceMesh>
    tensorSurfaceFieldSetType;

typedef FieldSetType<scalar, pointPatchField, pointMesh> scalarPointFieldSetType;
typedef FieldSetType<vector, pointPatchField, pointMesh> vectorPointFieldSetType;
typedef FieldSetType<symmTensor, pointPatchField, pointMesh>
    symmTensorPointFieldSetType;
typedef FieldSetType<sphericalTensor, pointPatchField, pointMesh>
    sphericalTensorPointFieldSetType;
typedef FieldSetType<tensor, pointPatchField, pointMesh> tensorPointFieldSetType;

}

#endif
// ************************************************************************* //

