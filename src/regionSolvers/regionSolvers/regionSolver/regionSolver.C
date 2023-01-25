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
#include "dynamicBlastFvMesh.H"

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
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    runTime_(mesh.time()),
    regions_(regions),
    dynMesh_(mesh),
    mesh_(dynMesh_),
    globalBoundary_(globalPolyBoundaryMesh::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolver::~regionSolver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam ::regionSolver::readControls
(
    const word& name,
    scalar& tol,
    scalar& relTol
) const
{
    if (!regions_.solutionControls().isDict(mesh_.name()))
    {
        return false;
    }
    const dictionary& regionDict =
        regions_.solutionControls().subDict(mesh_.name());
    if (regionDict.isDict(name))
    {
        regionDict.subDict(name).lookup("tolerance") >> tol;
        regionDict.subDict(name).lookup("relTol") >> relTol;
        return true;
    }
    return false;
}


void Foam::regionSolver::initialiseFields()
{}


void Foam::regionSolver::update()
{
    globalBoundary_.clearOut();
    globalBoundary_.update();
}


bool Foam::regionSolver::changeMesh()
{
    DebugInfo<< "Changing " << mesh_.name() << " mesh" << endl;
    return refineMesh(dynMesh_);
}


bool Foam::regionSolver::moveMesh(const IterType iter)
{
    DebugInfo<< "Moving " << mesh_.name() << " mesh" << endl;
    return dynMesh_.update();
}


void Foam::regionSolver::clear()
{}


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
