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

#include "NewmarkObjectMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionSolvers
{
    defineTypeNameAndDebug(Newmark, 0);
    addToRunTimeSelectionTable(objectMotionSolver, Newmark, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionSolvers::Newmark::Newmark
(
    const Newmark& solver,
    movingObject& body,
    objectMotionState& state,
    objectMotionState& state0
)
:
    objectMotionSolver(body, state, state0),
    gamma_(solver.gamma_),
    beta_(solver.beta_)
{}


Foam::objectMotionSolvers::Newmark::Newmark
(
    const dictionary& dict,
    movingObject& body,
    objectMotionState& state,
    objectMotionState& state0
)
:
    objectMotionSolver(body, state, state0),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 0.5)),
    beta_
    (
        max
        (
            0.25*sqr(gamma_ + 0.5),
            dict.lookupOrDefault<scalar>("beta", 0.25)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionSolvers::Newmark::~Newmark()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::objectMotionSolvers::Newmark::solve
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
    v() =
        tConstraints()
      & (v0() + deltaT*(gamma_*a() + (1 - gamma_)*a0()));

    // Correct angular momentum
    pi() =
        rConstraints()
      & (pi0() + deltaT*(gamma_*tau() + (1 - gamma_)*tau0()));

    // Correct position
    centreOfRotation() =
        centreOfRotation0()
      + (
            tConstraints()
          & (
                deltaT*v0()
              + sqr(deltaT)*(beta_*a() + (0.5 - beta_)*a0())
            )
        );

    // Correct orientation
    vector piDeltaT =
        rConstraints()
      & (
            deltaT*pi0()
          + sqr(deltaT)*(beta_*tau() + (0.5 - beta_)*tau0())
        );
    Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
    Q() = Qpi.first();
}


// ************************************************************************* //
