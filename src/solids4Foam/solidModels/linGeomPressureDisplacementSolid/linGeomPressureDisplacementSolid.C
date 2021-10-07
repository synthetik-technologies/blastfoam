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

#include "linGeomPressureDisplacementSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "linearElastic.H"
#include "findRefCell.H"
#include "adjustPhi.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomPressureDisplacementSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, linGeomPressureDisplacementSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomPressureDisplacementSolid::linGeomPressureDisplacementSolid
(
    dynamicFvMesh& mesh
)
:
    solidModel(typeName, mesh, nonLinGeom(), incremental()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    k_(mechanical().bulkModulus()),
    rK_(1.0/k_),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gradp_(fvc::grad(p_)),
    solvePressureEquationImplicitly_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "solvePressureEquationImplicitly",
            false
        )
    ),
    pressureRhieChowScaleFac_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "pressureRhieChowScaleFactor", 0.5
        )
    ),
    nInnerCorr_
    (
        solidModelDict().lookupOrDefault<int>("nInnerCorrectors", nCorr())
    )
{
    // Force p oldTime to be stored
    p_.oldTime();

    Info<< type() << ": solvePressureEquationImplicitly = "
        << solvePressureEquationImplicitly_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomPressureDisplacementSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        int iCorr = 0;
        SolverPerformance<vector> solverPerfD;
        SolverPerformance<vector>::debug = 0;

        Info<< "Solving the momentum equation for D and p" << endl;

        // Loop around displacement and pressure equations
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Momentum equation in terms of displacement
            // Note: sigma contains pressure term
            fvVectorMatrix DEqn
            (
                fvm::d2dt2(rho(), D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(D(), gradD(), impK_)
            );

            // Under-relaxation the linear system
            DEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DEqn);

            // Solve the linear system
            solverPerfD = DEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Store reciprocal of diagonal
            const surfaceScalarField rAUf(fvc::interpolate(1.0/DEqn.A()));

            // Calculate hydrostatic pressure
            if (solvePressureEquationImplicitly_)
            {
                int iInnerCorr = 0;
                SolverPerformance<scalar> solverPerfP;

                do
                {
                    // Store fields for under-relaxation and residual
                    // calculation
                    p_.storePrevIter();

                    // Pressure equation
                    fvScalarMatrix pEqn
                    (
                        fvm::Sp(rK_, p_)
                      - fvm::laplacian(rAUf, p_, "laplacian(Dp,p)")
                      + fvc::laplacian(rAUf, p_, "laplacian(Dp,p)")
                     ==
                      - fvc::div(D())
                      + pressureRhieChowScaleFac_
                       *(
                           fvc::laplacian(rAUf, p_, "laplacian(Dp,p)")
                         - fvc::div(rAUf*mesh().Sf() & fvc::interpolate(gradp_))
                        )
                    );

                    // Under-relaxation the linear system
                    pEqn.relax();

                    // Solve the linear system
                    solverPerfP = pEqn.solve();

                    // Under-relax the field
                    p_.relax();

                    // Update the gradient of pressure
                    gradp_ = fvc::grad(p_);

                    // Calculate the stress using run-time selectable mechanical
                    // law
                    mechanical().correct(sigma());

                    // Replace hydrostatic component of stress tensor
                    sigma() = dev(sigma()) - p_*I;
                }
                while
                (
                    !converged
                    (
                        iInnerCorr,
                        solverPerfP.initialResidual(),
                        solverPerfP.nIterations(),
                        p_,
                        false
                    )
                 && ++iInnerCorr < nInnerCorr_
                );
            }
            else
            {
                // Explicitly add Rhie-Chow smoothing to pressure field
                p_ = k_*
                (
                  - fvc::div(D())
                  + pressureRhieChowScaleFac_
                   *(
                       fvc::laplacian(rAUf, p_, "laplacian(Dp,p)")
                     - fvc::div(rAUf*mesh().Sf() & fvc::interpolate(gradp_))
                    )
                );

                // Update the gradient of pressure
                gradp_ = fvc::grad(p_);

                // Calculate the stress using run-time selectable mechanical
                // law
                mechanical().correct(sigma());

                // Replace hydrostatic component of stress tensor
                sigma() = dev(sigma()) - p_*I;
            }
        }
        while
        (
            !converged
            (
                iCorr,
                mag(solverPerfD.initialResidual()),
                max
                (
                    solverPerfD.nIterations()[0],
                    max
                    (
                        solverPerfD.nIterations()[1],
                        solverPerfD.nIterations()[2]
                    )
                ),
                D()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    return true;
}


tmp<vectorField> linGeomPressureDisplacementSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - n*pressure)
                  - (n & (pSigma - impK*pGradD))
                )*rImpK
            )
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
