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

#include "blastSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(blast, 0);
    addToRunTimeSelectionTable(solver, blast, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::blast::blast(fvMesh& mesh)
:
    explicitSolver(mesh),
    integrator_(mesh, false),
    fluidPtr_(compressibleSystem::New(mesh))
{
    integrator_.addSystem(fluidPtr_());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::blast::~blast()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::blast::CoNum() const
{
    return fluidPtr_->CoNum();
}


Foam::scalar Foam::solvers::blast::DiNum() const
{
    return fluidPtr_->DiNum();
}


void Foam::solvers::blast::solveExplicit()
{
    Info<< "Calculating Fluxes" << endl;
    integrator_.integrate(false); // No implicit
}


void Foam::solvers::blast::solveImplicit()
{
    integrator_.solveImplicit();
}


void Foam::solvers::blast::postSolve()
{
    Info<< "max(p): " << max(fluidPtr_->p()).value()
        << ", min(p): " << min(fluidPtr_->p()).value() << nl
        << "max(T): " << max(fluidPtr_->T()).value()
        << ", min(T): " << min(fluidPtr_->T()).value() << endl;

    integrator_.clear();
}


// ************************************************************************* //
