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

#include "blastParcelSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(blastParcel, 0);
    addToRunTimeSelectionTable(solver, blastParcel, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::blastParcel::blastParcel(fvMesh& mesh)
:
    explicitSolver(mesh),
    integrator_(mesh, false),
    fluidPtr_(coupledCompressibleSystem::New(mesh)),
    clouds_
    (
        parcelClouds::New
        (
            mesh_,
            fluidPtr_->rho(),
            fluidPtr_->U(),
            g_,
            fluidPtr_->thermo()
        )
    ),
    alphac_
    (
        IOobject
        (
            IOobject::groupName("alpha", parcelCloudList::cloudNamesName),
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        clouds_.alpha()
    )
{
    integrator_.addSystem(fluidPtr_());
    fluidPtr_->setDispersedVolumeFraction(alphac_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::blastParcel::~blastParcel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::blastParcel::CoNum() const
{
    return fluidPtr_->CoNum();
}


Foam::scalar Foam::solvers::blastParcel::DiNum() const
{
    return fluidPtr_->DiNum();
}


void Foam::solvers::blastParcel::preSolve()
{
    integrator_.preUpdateMesh();
    if (mesh_.dynamic() || mesh_.distributing())
    {
        clouds_.preUpdateMesh();
    }
    explicitSolver::preSolve();
}

void Foam::solvers::blastParcel::solveExplicit()
{
    fluidPtr_->decode();
    clouds_.evolve();
    alphac_ = clouds_.alpha();

    fluidPtr_->eSource() = clouds_.Sh(fluidPtr_->he());
    fluidPtr_->dragSource() = clouds_.SU(fluidPtr_->U());

    Info<< "Calculating Fluxes" << endl;
    integrator_.integrate(false);
}


void Foam::solvers::blastParcel::solveImplicit()
{
    integrator_.solveImplicit();
}


void Foam::solvers::blastParcel::postSolve()
{
    Info<< "max(p): " << max(fluidPtr_->p()).value()
        << ", min(p): " << min(fluidPtr_->p()).value() << nl
        << "max(T): " << max(fluidPtr_->T()).value()
        << ", min(T): " << min(fluidPtr_->T()).value() << endl;

    integrator_.clear();
}


// ************************************************************************* //
