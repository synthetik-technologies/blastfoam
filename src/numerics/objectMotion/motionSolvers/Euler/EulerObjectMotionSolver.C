/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "EulerObjectMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionSolvers
{
    defineTypeNameAndDebug(Euler, 0);
    addToRunTimeSelectionTable(objectMotionSolver, Euler, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionSolvers::Euler::Euler
(
    const Euler& solver,
    movingObject& body,
    objectMotionState& state,
    objectMotionState& state0
)
:
    objectMotionSolver(body, state, state0)
{}


Foam::objectMotionSolvers::Euler::Euler
(
    const dictionary& dict,
    movingObject& body,
    objectMotionState& state,
    objectMotionState& state0
)
:
    objectMotionSolver(body, state, state0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionSolvers::Euler::~Euler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::objectMotionSolvers::Euler::solve
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    // Update the linear acceleration and torque
    updateAcceleration(fGlobal, tauGlobal);

    // Correct linear velocity
    v() = tConstraints() & (v0() + deltaT*a());

    // Correct angular momentum
    pi() = rConstraints() & (pi0() + deltaT*tau());

    // Correct position
    centreOfRotation() = centreOfRotation0() + deltaT*v();

    // Correct orientation
    Tuple2<tensor, vector> Qpi = rotate(Q0(), pi(), deltaT);
    Q() = Qpi.first();
}


// ************************************************************************* //
