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

#include "coupledSolidTractionFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"

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

Foam::regionSolvers::solid::solid
(
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    regionSolver(mesh, regions),
    solid_(solidModel::New(dynMesh_))
{
    // Add Displacement field to track error
    accelerationSchemes_.addField(solid_->D());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::solid::~solid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::solid::initialiseMesh(const IterType)
{}


void Foam::regionSolvers::solid::initialiseFields()
{
    solidTractionFvPatchVectorField::canRelax = false;

    const_cast<volVectorField&>
    (
        solid_->solutionD()
    ).correctBoundaryConditions();
    solid_->update();
}


void Foam::regionSolvers::solid::initialise()
{
    solidTractionFvPatchVectorField::canRelax = false;

    solid_->initialize();

    const_cast<volVectorField&>
    (
        solid_->solutionD()
    ).correctBoundaryConditions();
    solid_->update();
}


bool Foam::regionSolvers::solid::moveMesh(const IterType iter)
{
    if (iter == FINAL_ITER)
    {
        solidTractionFvPatchVectorField::canRelax = false;
    }
    else
    {
        solidTractionFvPatchVectorField::canRelax = true;
    }

    regionSolver::moveMesh(iter);
    return max(mag(solid_->DD())).value() > small;
}


void Foam::regionSolvers::solid::solve()
{
    SolverPerformance<vector>::debug = 0;

    solid_->D().storePrevIter();
    solid_->evolve();

    accelerationSchemes_.updateError();

    const volVectorField& D = solid_->solutionD();

    vector forceSum = Zero;
    Info<< "External forces:" << incrIndent << endl;
    forAll(D.boundaryField(), patchi)
    {
        const fvPatchVectorField& pD = D.boundaryField()[patchi];
        if (isA<coupledSolidTractionFvPatchVectorField>(pD))
        {
            const coupledSolidTractionFvPatchVectorField& cst =
                dynamicCast<const coupledSolidTractionFvPatchVectorField>(pD);
            forceSum += cst.force();
            Info<< indent << pD.patch().name() << ":" << nl << incrIndent
                << indent << "solid = " << cst.force() << nl
                << indent << "fluid = " << cst.forceNbr() << decrIndent << endl;
        }
        else if (isA<solidTractionFvPatchVectorField>(pD))
        {
            const solidTractionFvPatchVectorField& st =
                dynamicCast<const solidTractionFvPatchVectorField>(pD);
            forceSum += st.force();
            Info<< indent << pD.patch().name() << ": "
                << st.force() << endl;
        }
    }
    Info<< indent << "Total: " << forceSum << decrIndent << nl << endl;

    // Turn solver information back on
    SolverPerformance<vector>::debug = 1;
}


void Foam::regionSolvers::solid::clear(const bool full)
{
    regionSolver::clear(full);
    solidTractionFvPatchVectorField::canRelax = true;
    solid_->updateTotalFields();
}


Foam::scalar Foam::regionSolvers::solid::CoNum() const
{
    scalar CoNum = solid_->CoNum();
    Info<< mesh_.name() << ": Courant Number max: " << CoNum << endl;
    return CoNum;
}


Foam::scalar Foam::regionSolvers::solid::maxCo() const
{
    return
        min
        (
            runTime_.controlDict().lookupOrDefault
            (
                mesh_.name() + "MaxCo",
                runTime_.controlDict().lookup<scalar>("maxCo")
            ),
            solid_->maxCoNum()
        );
}


Foam::scalar Foam::regionSolvers::solid::newDeltaT() const
{
    return min(regionSolver::newDeltaT(), solid_->mechanical().newDeltaT());
}


// ************************************************************************* //
