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

#include "blastReactingSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(blastReacting, 0);
    addToRunTimeSelectionTable(solver, blastReacting, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::blastReacting::blastReacting(fvMesh& mesh)
:
    explicitSolver(mesh),
    integrator_(mesh),
    fluid_
    (
        reactingCompressibleSystem
        (
            IOdictionary
            (
                IOobject
                (
                    ::Foam::physicalProperties::typeName,
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            ),
            mesh
        )
    )
{
    integrator_.addSystem(fluid_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::blastReacting::~blastReacting()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::blastReacting::CoNum() const
{
    return fluid_.CoNum();
}


void Foam::solvers::blastReacting::solve()
{
    Info<< "Calculating Fluxes" << endl;
    integrator_.integrate();
}


void Foam::solvers::blastReacting::postSolve()
{
    Info<< "max(p): " << max(fluid_.p()).value()
        << ", min(p): " << min(fluid_.p()).value() << nl
        << "max(T): " << max(fluid_.T()).value()
        << ", min(T): " << min(fluid_.T()).value() << endl;

    integrator_.clear();
}


// ************************************************************************* //
