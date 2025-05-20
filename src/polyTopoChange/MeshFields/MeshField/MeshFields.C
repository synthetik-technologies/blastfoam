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

\*---------------------------------------------------------------------------*/

#include "MeshFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// defineTemplateTypeNameAndDebug(cellBoolMeshField, 0);
// defineTemplateTypeNameAndDebug(cellLabelMeshField, 0);
defineTemplateTypeNameAndDebug(cellScalarMeshField, 0);
defineTemplateTypeNameAndDebug(cellVectorMeshField, 0);
defineTemplateTypeNameAndDebug(cellSphericalTensorMeshField, 0);
defineTemplateTypeNameAndDebug(cellSymmTensorMeshField, 0);
defineTemplateTypeNameAndDebug(cellTensorMeshField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// defineTemplateTypeNameAndDebug(faceBoolMeshField, 0);
// defineTemplateTypeNameAndDebug(faceLabelMeshField, 0);
defineTemplateTypeNameAndDebug(faceScalarMeshField, 0);
defineTemplateTypeNameAndDebug(faceVectorMeshField, 0);
defineTemplateTypeNameAndDebug(faceSphericalTensorMeshField, 0);
defineTemplateTypeNameAndDebug(faceSymmTensorMeshField, 0);
defineTemplateTypeNameAndDebug(faceTensorMeshField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// defineTemplateTypeNameAndDebug(pointBoolMeshField, 0);
// defineTemplateTypeNameAndDebug(pointLabelMeshField, 0);
defineTemplateTypeNameAndDebug(pointScalarMeshField, 0);
defineTemplateTypeNameAndDebug(pointVectorMeshField, 0);
defineTemplateTypeNameAndDebug(pointSphericalTensorMeshField, 0);
defineTemplateTypeNameAndDebug(pointSymmTensorMeshField, 0);
defineTemplateTypeNameAndDebug(pointTensorMeshField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //