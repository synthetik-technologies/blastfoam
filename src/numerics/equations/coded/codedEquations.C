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

#include "CodedEquation.H"
#include "CodedUnivariateEquation.H"
#include "CodedMultivariateEquation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(CodedEquation<scalar>, 0);
addTemplatedToRunTimeSelectionTable
(
    equation,
    CodedEquation,
    scalar,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedEquation<vector>, 0);
addTemplatedToRunTimeSelectionTable
(
    equation,
    CodedEquation,
    vector,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedEquation<symmTensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    equation,
    CodedEquation,
    symmTensor,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedEquation<sphericalTensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    equation,
    CodedEquation,
    sphericalTensor,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedEquation<tensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    equation,
    CodedEquation,
    tensor,
    dictionary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(CodedUnivariateEquation<scalar>, 0);
addTemplatedToRunTimeSelectionTable
(
    univariateEquation,
    CodedUnivariateEquation,
    scalar,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedUnivariateEquation<vector>, 0);
addTemplatedToRunTimeSelectionTable
(
    univariateEquation,
    CodedUnivariateEquation,
    vector,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedUnivariateEquation<symmTensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    univariateEquation,
    CodedUnivariateEquation,
    symmTensor,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedUnivariateEquation<sphericalTensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    univariateEquation,
    CodedUnivariateEquation,
    sphericalTensor,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedUnivariateEquation<tensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    univariateEquation,
    CodedUnivariateEquation,
    tensor,
    dictionary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(CodedMultivariateEquation<scalar>, 0);
addTemplatedToRunTimeSelectionTable
(
    multivariateEquation,
    CodedMultivariateEquation,
    scalar,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedMultivariateEquation<vector>, 0);
addTemplatedToRunTimeSelectionTable
(
    multivariateEquation,
    CodedMultivariateEquation,
    vector,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedMultivariateEquation<symmTensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    multivariateEquation,
    CodedMultivariateEquation,
    symmTensor,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedMultivariateEquation<sphericalTensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    multivariateEquation,
    CodedMultivariateEquation,
    sphericalTensor,
    dictionary
);

defineTemplateTypeNameAndDebug(CodedMultivariateEquation<tensor>, 0);
addTemplatedToRunTimeSelectionTable
(
    multivariateEquation,
    CodedMultivariateEquation,
    tensor,
    dictionary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
