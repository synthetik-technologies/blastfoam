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

#include "unsNonLinGeomTotalLagSolid.H"
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

defineTypeNameAndDebug(unsNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, unsNonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


scalar unsNonLinGeomTotalLagSolid::residual(const volVectorField& D) const
{
    return
        gMax
        (
            DimensionedField<double, volMesh>
            (
                mag(D.internalField() - D.prevIter().internalField())
               /max
                (
                    gMax
                    (
                        DimensionedField<double, volMesh>
                        (
                            mag
                            (
                                D.internalField() - D.oldTime().internalField()
                            )
                        )
                    ),
                    SMALL
                )
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsNonLinGeomTotalLagSolid::unsNonLinGeomTotalLagSolid(dynamicFvMesh& mesh)
:
    unsTotalLagSolid<unsTotalDispSolid>(typeName, mesh),
    K_
    (
        solidModelDict().lookupOrDefault<dimensionedScalar>
        (
            "K",
            dimensionedScalar("K", dimless/dimTime, 0)
        )
    ),
    relativeTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "solutionTolerance",
            solutionTol()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool unsNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    scalar initialResidual = 0;
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<vector>::debug = 0;
    scalar res = 1.0;
    scalar maxRes = 0;
    scalar curConvergenceTolerance = solutionTol();

    // Reset enforceLinear switch
    enforceLinear() = false;

    do
    {
        if (SolverPerformance<vector>::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        // Store previous iteration to allow under-relaxation and residual
        // calculation
        D().storePrevIter();

        // Construct momentum equation in total Lagrangian form where gradients
        // are calculated directly at the faces
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
          + rho()*g()
        );

        // Add damping
        if (K_.value() > SMALL)
        {
            DEqn += K_*rho()*fvm::ddt(D());
        }

        // Enforce linear to improve convergence
        if (enforceLinear())
        {
            // Replace nonlinear terms with linear
            // Note: the mechanical law could still be nonlinear
            DEqn +=
                fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
              - fvc::div(mesh().Sf() & sigmaf_);
        }

        // Under-relax the linear system
        DEqn.relax();

        // Solve the system
        solverPerfD = DEqn.solve();

        // Under-relax displacement field
        relaxField(D(), iCorr);

        if (iCorr == 0)
        {
            initialResidual = mag(solverPerfD.initialResidual());
        }

        update();

        // Check if outer loops are diverging
        if (!enforceLinear())
        {
            checkEnforceLinear(Jf_);
        }

        // Calculate relative momentum residual
        res = residual(D());

        if (res > maxRes)
        {
            maxRes = res;
        }

        curConvergenceTolerance = maxRes*relativeTol_;

        if (curConvergenceTolerance < solutionTol())
        {
            curConvergenceTolerance = solutionTol();
        }

        if
        (
            SolverPerformance<vector>::debug
         || (iCorr % infoFrequency()) == 0
         || res < curConvergenceTolerance
         || maxIterReached() == nCorr()
        )
        {
            Info<< "Corr " << iCorr << ", relative residual = " << res << endl;
        }

        if (maxIterReached() == nCorr())
        {
            maxIterReached()++;
        }
    }
    while (res > curConvergenceTolerance && ++iCorr < nCorr());

    // Velocity
    U() = fvc::ddt(D());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Print summary of residuals
    Info<< solverPerfD.solverName() << ": Solving for " << D().name()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerfD.initialResidual()
        << ", No outer iterations = " << iCorr << nl
        << " Max relative residual = " << maxRes
        << ", Relative residual = " << res
        << ", enforceLinear = " << enforceLinear() << endl;

    SolverPerformance<vector>::debug = 1;

    if (enforceLinear())
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
