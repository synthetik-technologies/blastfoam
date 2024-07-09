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

#include "blastRegionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(blast, 0);
    addToRunTimeSelectionTable(regionSolver, blast, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::blast::blast
(
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    fluid(mesh, regions),
    g_
    (
        IOobject
        (
            "g",
            runTime_.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
    integrator_(mesh_),
    fluid_(compressibleSystem::New(mesh_))
{
    integrator_.addSystem(fluid_());
}


Foam::regionSolvers::blast::blast
(
    dynamicFvMesh& mesh,
    const word& type,
    const regionSolverList& regions
)
:
    fluid(mesh, regions),
    g_
    (
        IOobject
        (
            "g",
            runTime_.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
    integrator_(mesh_),
    fluid_(compressibleSystem::New(type, mesh_))
{
    integrator_.addSystem(fluid_());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::blast::~blast()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionSolvers::blast::changeMesh()
{
    integrator_.preUpdateMesh();
    return fluid::changeMesh();
}


void Foam::regionSolvers::blast::solve()
{
    integrator_.integrate();

    Info<< "max(p): " << max(fluid_->p()).value()
        << ", min(p): " << min(fluid_->p()).value() << endl;
    Info<< "max(T): " << max(fluid_->T()).value()
        << ", min(T): " << min(fluid_->T()).value() << endl;

    integrator_.clear();
}


Foam::scalar Foam::regionSolvers::blast::CoNum() const
{
    return fluid_->CoNum();
}


Foam::scalar Foam::regionSolvers::blast::maxCo() const
{
    return min(fluid::maxCo(), integrator_.maxCo());
}

// ************************************************************************* //
