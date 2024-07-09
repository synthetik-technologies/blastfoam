/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "polyMesh.H"
#include "feFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplate2TypeNameAndDebug(feScalarField::Internal, 0);
defineTemplate2TypeNameAndDebug(feVectorField::Internal, 0);
defineTemplate2TypeNameAndDebug
(
    feSphericalTensorField::Internal,
    0
);
defineTemplate2TypeNameAndDebug
(
    feSymmTensorField::Internal,
    0
);
defineTemplate2TypeNameAndDebug(feTensorField::Internal, 0);


defineTemplateTypeNameAndDebug(feScalarField, 0);
defineTemplateTypeNameAndDebug(feVectorField, 0);
defineTemplateTypeNameAndDebug(feSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(feSymmTensorField, 0);
defineTemplateTypeNameAndDebug(feTensorField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
