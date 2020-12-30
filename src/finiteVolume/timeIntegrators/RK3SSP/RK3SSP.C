/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "RK3SSP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK3SSP, 0);
    addToRunTimeSelectionTable(timeIntegrator, RK3SSP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK3SSP::RK3SSP
(
    phaseCompressibleSystem& fluid
)
:
    timeIntegrator(fluid)
{
    fluid_.setODEFields
    (
        3,
        {true, false, false},
        {false, false, false}
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK3SSP::~RK3SSP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RK3SSP::integrate()
{
    // Update and store original fields
    fluid_.update();
    fluid_.solve(1, {1.0}, {1.0});

    // Update and store 1st step
    fluid_.update();
    fluid_.solve(2, {0.75, 0.25}, {0.0, 0.25});

    // Update and store 1st step
    fluid_.update();
    fluid_.solve(3, {1.0/3.0, 0.0, 2.0/3.0}, {0.0, 0.0, 2.0/3.0});
}
// ************************************************************************* //
