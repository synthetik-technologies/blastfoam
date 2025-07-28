/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

InClass
    Foam::MeshField

Description

SourceFiles
    MeshField.C

\*---------------------------------------------------------------------------*/

#ifndef MeshFieldsFwd_H
#define MeshFieldsFwd_H

#include "MeshField.H"
#include "MeshFieldUpdatersFwd.H"
#include "fieldTypes.H"
#include "cellGeoMesh.H"
#include "faceGeoMesh.H"
#include "pointGeoMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// typedef MeshField<bool, cellGeoMesh> cellBoolMeshField;
// typedef MeshField<label, cellGeoMesh> cellLabelMeshField;
typedef MeshField<scalar, cellGeoMesh> cellScalarMeshField;
typedef MeshField<vector, cellGeoMesh> cellVectorMeshField;
typedef MeshField<sphericalTensor, cellGeoMesh> cellSphericalTensorMeshField;
typedef MeshField<symmTensor, cellGeoMesh> cellSymmTensorMeshField;
typedef MeshField<tensor, cellGeoMesh> cellTensorMeshField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// typedef MeshField<bool, faceGeoMesh> faceBoolMeshField;
// typedef MeshField<label, faceGeoMesh> faceLabelMeshField;
typedef MeshField<scalar, faceGeoMesh> faceScalarMeshField;
typedef MeshField<vector, faceGeoMesh> faceVectorMeshField;
typedef MeshField<sphericalTensor, faceGeoMesh> faceSphericalTensorMeshField;
typedef MeshField<symmTensor, faceGeoMesh> faceSymmTensorMeshField;
typedef MeshField<tensor, faceGeoMesh> faceTensorMeshField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// typedef MeshField<bool, pointGeoMesh> pointBoolMeshField;
// typedef MeshField<label, pointGeoMesh> pointLabelMeshField;
typedef MeshField<scalar, pointGeoMesh> pointScalarMeshField;
typedef MeshField<vector, pointGeoMesh> pointVectorMeshField;
typedef MeshField<sphericalTensor, pointGeoMesh> pointSphericalTensorMeshField;
typedef MeshField<symmTensor, pointGeoMesh> pointSymmTensorMeshField;
typedef MeshField<tensor, pointGeoMesh> pointTensorMeshField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //