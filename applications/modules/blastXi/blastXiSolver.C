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

#include "blastXiSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(blastXi, 0);
    addToRunTimeSelectionTable(solver, blastXi, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::blastXi::blastXi(fvMesh& mesh)
:
    explicitSolver(mesh),
    integrator_(mesh, false),
    fluid_
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
{
    integrator_.addSystem(fluid_);
    fluxSchemeBase::needEnergyFlux = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::blastXi::~blastXi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::blastXi::CoNum() const
{
    return fluid_.CoNum();
}



Foam::scalar Foam::solvers::blastXi::DiNum() const
{
    return fluid_.DiNum();
}


void Foam::solvers::blastXi::solveExplicit()
{
    Info<< "Calculating Fluxes" << endl;
    integrator_.integrate
    (
        true,   // doExplicit
        true,   // doStore
        false,  // doImplicit
        false,  // doPost
        false   // doClear
    );
}

void Foam::solvers::blastXi::solveImplicit()
{
    integrator_.solveImplicit();
}


void Foam::solvers::blastXi::postSolve()
{
    Info<< "max(p): " << max(fluid_.p()).value()
        << ", min(p): " << min(fluid_.p()).value() << nl
        << "max(T): " << max(fluid_.T()).value()
        << ", min(T): " << min(fluid_.T()).value() << endl;

    integrator_.postUpdate();
    integrator_.clear();
}


// ************************************************************************* //
