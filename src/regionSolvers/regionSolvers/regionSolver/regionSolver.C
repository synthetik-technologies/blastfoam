/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "regionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSolver, 0);
    defineRunTimeSelectionTable(regionSolver, dictionary);
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolver::regionSolver
(
    fvMesh& mesh,
    const regionSolverList& regions
)
:
    runTime_(mesh.time()),
    regions_(regions),
    mesh_(mesh),
    globalBoundary_(globalPolyBoundaryMesh::New(mesh)),
    accelerationSchemes_
    (
        mesh_,
        regions_.regionProperties(),
        "solutionControls"
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolver::~regionSolver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolver::storePrevIter()
{
    accelerationSchemes_.storePrevIter();
}


void Foam::regionSolver::initialiseFields()
{
    // Look up fields to relax
}


void Foam::regionSolver::update()
{
    globalBoundary_.clearOut();
    globalBoundary_.update();
}


bool Foam::regionSolver::changeMesh()
{
    DebugInfo<< "Changing " << mesh_.name() << " mesh" << endl;
    if (mesh_.update())
    {
        this->clear(true);
        return true;
    }
    return false;
}


bool Foam::regionSolver::moveMesh(const IterType iter)
{
    DebugInfo<< "Moving " << mesh_.name() << " mesh" << endl;
    return mesh_.move();
}


void Foam::regionSolver::clear(const bool full)
{
    DebugInfo<< "Clearing " << mesh_.name() << endl;
    accelerationSchemes_.clear(full);
    globalBoundary_.write();
}


Foam::scalar Foam::regionSolver::maxCo() const
{
    return runTime_.controlDict().lookupOrDefault
        (
            mesh_.name() + "MaxCo",
            runTime_.controlDict().lookup<scalar>("maxCo")
        );
}

Foam::scalar Foam::regionSolver::newDeltaT() const
{
    scalar maxDeltaTFact =
        this->maxCo()/(this->CoNum() + small);
    scalar deltaTFact =
        min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    return deltaTFact*runTime_.deltaTValue();
}

// ************************************************************************* //
