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

#include "RK4SSPTimeIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK4SSP, 0);
    addToRunTimeSelectionTable(timeIntegrator, RK4SSP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK4SSP::RK4SSP
(
    const fvMesh& mesh
)
:
    timeIntegrator(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK4SSP::~RK4SSP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RK4SSP::setODEFields(integrationSystem& system)
{
    system.setODEFields
    (
        4,
        {true, true, true, false},
        {true, true, true, false}
    );
}


void Foam::timeIntegrators::RK4SSP::integrate()
{
    // Update and store original fields
    forAll(systems_, i)
    {
        systems_[i].update();
        systems_[i].solve(1, {1.0}, {0.5});
    }

    // Update and store 1st step
    forAll(systems_, i)
    {
        systems_[i].update();
        systems_[i].solve(2, {a20, a21}, {b20, b21});
    }

    // Update and store 1st step
    forAll(systems_, i)
    {
        systems_[i].update();
        systems_[i].solve(3, {a30, a31, a32}, {b30, b31, b32});
    }

    // Update and store 1st step
    forAll(systems_, i)
    {
        systems_[i].update();
        systems_[i].solve
        (
            4,
            {a40, a41, a42, a43},
            {b40, b41, b42, b43}
        );
    }
}
// ************************************************************************* //
