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

#include "solidRegionSolver.H"
#include "globalPolyBoundaryMesh.H"
#include "SolverPerformance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(solid, 0);
    addToRunTimeSelectionTable(regionSolver, solid, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::solid::solid(dynamicFvMesh& mesh)
:
    regionSolver(mesh),
    solid_(solidModel::New(dynMesh_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::solid::~solid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Initialise the mesh
void Foam::regionSolvers::solid::initialiseMesh()
{}

//- Initialise the solver
void Foam::regionSolvers::solid::initialise()
{
    solid_->initialize();
}

//- Solve the model
void Foam::regionSolvers::solid::solve()
{
    SolverPerformance<vector>::debug = 0;

    solid_->evolve();
    solid_->updateTotalFields();

    // Turn solver information back on
    SolverPerformance<vector>::debug = 1;

    //- Clear global Patches since displacement may have changed
    globalPolyBoundaryMesh::New(mesh_).movePoints();
}

//- Return the Courant number
Foam::scalar Foam::regionSolvers::solid::CoNum() const
{
    return solid_->CoNum();
}

//- Return the maximum Courant number
Foam::scalar Foam::regionSolvers::solid::maxCo() const
{
    return solid_->maxCoNum();
}

// ************************************************************************* //
