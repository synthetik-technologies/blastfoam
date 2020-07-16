/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "CarnahanStarlingRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(CarnahanStarling, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        CarnahanStarling,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::CarnahanStarling::CarnahanStarling
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    radialModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::CarnahanStarling::~CarnahanStarling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::CarnahanStarling::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{

    return
        1.0/(1 - phase1)
        + 3*phase1/(2*sqr(1 - phase1))
        + sqr(phase1)/(2*pow3(1 - phase1));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::CarnahanStarling::gs0prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        2.5/sqr(1 - phase1)
        + 4*phase1/pow3(1 - phase1)
        + 1.5*sqr(phase1)/pow4(1 - phase1);
}


// ************************************************************************* //
