/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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

#include "explicitSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(explicitSolver, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::explicitSolver::explicitSolver(fvMesh& mesh)
:
    solver(mesh),
    maxCo_(0.5),
    maxDeltaT_(vGreat),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector::zero)
    )
{
    steady = false;
    LTS = false;

    const_cast<dictionary&>(pimple.dict()).set("nOuterCorrectors", 1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::explicitSolver::~explicitSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solvers::explicitSolver::read()
{
    const_cast<dictionary&>(pimple.dict()).set("nOuterCorrectors", 1);

    solver::read();

    maxCo_ =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.5);

    maxDeltaT_ =
        runTime.controlDict().found("maxDeltaT")
      ? runTime.controlDict().lookup<scalar>("maxDeltaT", runTime.userUnits())
      : vGreat;

    return true;
}


Foam::scalar Foam::solvers::explicitSolver::maxDeltaT() const
{
    const scalar Co = this->CoNum();

    scalar deltaT = min(fvModels().maxDeltaT(), maxDeltaT_);

    if (maxCo_ < vGreat && Co > small)
    {
        deltaT = min(deltaT, maxCo_/Co*runTime.deltaTValue());
    }

    return deltaT;
}


void Foam::solvers::explicitSolver::preSolve()
{
    mesh_.update();
}


void Foam::solvers::explicitSolver::moveMesh()
{
    if (pimple.firstIter())
    {
        mesh_.move();
    }
}


void Foam::solvers::explicitSolver::motionCorrector()
{}


void Foam::solvers::explicitSolver::prePredictor()
{}

void Foam::solvers::explicitSolver::momentumPredictor()
{
    Info<<runTime.deltaTValue()<<endl;
    this->solve();
}

void Foam::solvers::explicitSolver::thermophysicalPredictor()
{}


void Foam::solvers::explicitSolver::pressureCorrector()
{}


void Foam::solvers::explicitSolver::postCorrector()
{}


void Foam::solvers::explicitSolver::postSolve()
{}

// ************************************************************************* //
