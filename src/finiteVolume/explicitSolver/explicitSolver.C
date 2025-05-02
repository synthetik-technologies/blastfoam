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
    defineTypeNameAndDebug(explicitSolver, 0);\
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitSolver::explicitSolver(fvMesh& mesh)
:
    solver(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::explicitSolver::~explicitSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::explicitSolver::preSolve()
{
    mesh_.update();
}


void Foam::explicitSolver::moveMesh()
{
    if (pimple.firstIter())
    {
        mesh_.move();
    }
}


void Foam::explicitSolver::motionCorrector()
{}


void Foam::explicitSolver::prePredictor()
{}

void Foam::explicitSolver::momentumPredictor()
{
    this->solve();
}

void Foam::explicitSolver::thermophysicalPredictor()
{}


void Foam::explicitSolver::pressureCorrector()
{}


void Foam::explicitSolver::postCorrector()
{}


void Foam::explicitSolver::postSolve()
{}


// ************************************************************************* //
