/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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

#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
defineNamedTemplateTypeNameAndDebug(ReconstructionScheme<scalar>, 0);
defineTemplateRunTimeSelectionTable
(
    ReconstructionScheme<scalar>,
    dictionary
);

defineNamedTemplateTypeNameAndDebug(ReconstructionScheme<vector>, 0);
defineTemplateRunTimeSelectionTable
(
    ReconstructionScheme<vector>,
    dictionary
)

defineNamedTemplateTypeNameAndDebug(ReconstructionScheme<symmTensor>, 0);
defineTemplateRunTimeSelectionTable
(
    ReconstructionScheme<symmTensor>,
    dictionary
)

defineNamedTemplateTypeNameAndDebug(ReconstructionScheme<sphericalTensor>, 0);
defineTemplateRunTimeSelectionTable
(
    ReconstructionScheme<sphericalTensor>,
    dictionary
)

defineNamedTemplateTypeNameAndDebug(ReconstructionScheme<tensor>, 0);
defineTemplateRunTimeSelectionTable
(
    ReconstructionScheme<tensor>,
    dictionary
)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
