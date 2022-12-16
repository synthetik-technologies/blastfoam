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

#include "hydrostaticRegionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(hydrostatic, 0);
    addToRunTimeSelectionTable(regionSolver, hydrostatic, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::hydrostatic::hydrostatic
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
    fluid_(compressibleSystem::New(mesh_)),
    atmosphereProperties_
    (
        IOobject
        (
            "atmosphereProperties",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    atmosphere_(atmosphereModel::New(mesh, atmosphereProperties_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::hydrostatic::~hydrostatic()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::hydrostatic::solve()
{
    Info<< "Setting hydrostatic pressure" << endl;
    atmosphere_->createAtmosphere
    (
        dynamicCast<fluidBlastThermo>(fluid_->thermo())
    );

    Info<< "max(p): " << max(fluid_->p()).value()
        << ", min(p): " << min(fluid_->p()).value() << endl;
    Info<< "max(T): " << max(fluid_->T()).value()
        << ", min(T): " << min(fluid_->T()).value() << endl;
}


Foam::scalar Foam::regionSolvers::hydrostatic::CoNum() const
{
    return 0.0;
}


Foam::scalar Foam::regionSolvers::hydrostatic::maxCo() const
{
    return great;
}

// ************************************************************************* //
