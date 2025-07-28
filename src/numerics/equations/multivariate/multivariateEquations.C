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

#include "MultivariateEquationsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(multivariateEquation<scalar>, 0);
defineTemplateRunTimeSelectionTable(multivariateEquation<scalar>, dictionary);

defineTemplateTypeNameAndDebug(multivariateEquation<vector>, 0);
defineTemplateRunTimeSelectionTable(multivariateEquation<vector>, dictionary);

defineTemplateTypeNameAndDebug(multivariateEquation<symmTensor>, 0);
defineTemplateRunTimeSelectionTable(multivariateEquation<symmTensor>, dictionary);

defineTemplateTypeNameAndDebug(multivariateEquation<sphericalTensor>, 0);
defineTemplateRunTimeSelectionTable(multivariateEquation<sphericalTensor>, dictionary);

defineTemplateTypeNameAndDebug(multivariateEquation<tensor>, 0);
defineTemplateRunTimeSelectionTable(multivariateEquation<tensor>, dictionary);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
