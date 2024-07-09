/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "unsUpdatedLagrangianSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsUpdatedLagrangianSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, unsUpdatedLagrangianSolid, dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsUpdatedLagrangianSolid::unsUpdatedLagrangianSolid(dynamicFvMesh& mesh)
:
    UnsUpdatedLagrangianGeomSolid<unsIncrementalSolid>(typeName, mesh)
{
    //- Dummy Call to make sure the necessary old fields are initialized
    fvc::d2dt2(rho().oldTime(), D().oldTime());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool unsUpdatedLagrangianSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    SolverPerformance<vector> solverPerfDD;
    SolverPerformance<vector>::debug = 0;

    Info<< "Solving the momentum equation for DD" << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Momentum equation incremental updated Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho(), DD())
          + fvc::d2dt2(rho().oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")

          //- Relative Piola (relative to current configuration)
          + fvc::div(this->Pf() & mesh().Sf())

          + rho()*g()
        );

        // Under-relax the linear system
        DDEqn.relax();

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field
        relaxField(DD(), iCorr);

        update();
    }
    while
    (
       !converged
        (
            iCorr,
            mag(solverPerfDD.initialResidual()),
            max
            (
                solverPerfDD.nIterations()[0],
                max
                (
                    solverPerfDD.nIterations()[1],
                    solverPerfDD.nIterations()[2]
                )
            ),
            DD()
        ) && ++iCorr < nCorr()
    );

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
