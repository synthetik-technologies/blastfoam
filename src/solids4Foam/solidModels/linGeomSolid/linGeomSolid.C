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

#include "linGeomSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomSolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomSolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomSolid::linGeomSolid(dynamicFvMesh& mesh)
:
    solidModel(typeName, mesh, nonLinGeom(), incremental()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{
    DDisRequired();

    //- Dummy Call to make sure the necessary old fields are initialized
    fvc::d2dt2(rho().oldTime(), D().oldTime());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        int iCorr = 0;
        SolverPerformance<vector> solverPerfDD;
        SolverPerformance<vector>::debug = 0;

        Info<< "Solving the momentum equation for DD" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            DD().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DDEqn
            (
                rho()*fvm::d2dt2(DD())
              + rho()*fvc::d2dt2(D().oldTime())
             == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
              - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + mechanical().RhieChowCorrection(DD(), gradDD())
            );

            // Under-relaxation the linear system
            DDEqn.relax();

            // Solve the linear system
            solverPerfDD = DDEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(DD(), iCorr);

            // Update the total displacement
            D() = D().oldTime() + DD();

            // Update gradient of displacement increment
            mechanical().grad(DD(), gradDD());

            // Update gradient of total displacement
            gradD() = gradD().oldTime() + gradDD();

            // Calculate the stress using run-time selectable mechanical law
            const volScalarField DDEqnA("DDEqnA", DDEqn.A());
            mechanical().correct(sigma());
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

        // Update point displacement increment
        mechanical().interpolate(DD(), pointDD());

        // Update point displacement
        pointD() = pointD().oldTime() + pointDD();

        // Update velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    return true;
}


tmp<vectorField> linGeomSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& pImpK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& pRImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField pN(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - pN*pressure)
              - (pN & (pSigma - pImpK*pGradDD))
            )*pRImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
