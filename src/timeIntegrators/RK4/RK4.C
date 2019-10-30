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

#include "RK4.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK4, 0);
    addToRunTimeSelectionTable(timeIntegrator, RK4, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK4::RK4
(
    phaseCompressibleSystem& fluid
)
:
    timeIntegrator(fluid)
{
    fluid_.setODEFields
    (
        4,
        {true, false, false, false},
        {true, true, true, false}
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK4::~RK4()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RK4::integrate()
{
    // Update and store original fields
    fluid_.update();
    fluid_.solve(1, {1.0}, {0.5});

    // Update and store 1st step
    fluid_.update();
    fluid_.solve(2, {0.0, 1.0}, {0.0, 0.5});

    // Update and store 1st step
    fluid_.update();
    fluid_.solve(3, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0});

    // Update and store 1st step
    fluid_.update();
    fluid_.solve
    (
        4,
        {1.0, 0.0, 0.0, 0.0},
        {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}
    );
}
// ************************************************************************* //
