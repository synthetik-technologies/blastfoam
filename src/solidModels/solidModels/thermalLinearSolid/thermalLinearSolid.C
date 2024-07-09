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

#include "thermalLinearSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalLinearSolid, 0);
addToRunTimeSelectionTable(solidModel, thermalLinearSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool thermalLinearSolid::converged
(
    const int iCorr,
    const SolverPerformance<vector>& solverPerfD,
    const SolverPerformance<scalar>& solverPerfHE,
    const volVectorField& D,
    const volScalarField& he
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar absResidualHE =
        gMax
        (
            DimensionedField<double, volMesh>
            (
                mag(he.internalField() - he.prevIter().internalField())
            )
        );
    const scalar residualHE =
        absResidualHE
       /max
        (
            gMax
            (
                DimensionedField<double, volMesh>
                (
                    mag(he.internalField() - he.oldTime().internalField())
                )
            ),
            SMALL
        );

    const scalar residualD =
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
                            mag(D.internalField() - D.oldTime().internalField())
                        )
                    ),
                    SMALL
                )
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();
    const scalar materialRelResidual = this->mechanical().relResidual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if
    (
        iCorr > 1
     && (
            materialResidual < materialTol()
         || materialRelResidual < materialRelTol()
        )
    )
    {
        bool convergedD = false;
        bool convergedHE = false;

        if
        (
            (
                mag(solverPerfD.initialResidual()) < solutionTol()
             && residualD < solutionTol()
            )
         || mag(solverPerfD.initialResidual()) < alternativeTol()
         || residualD < alternativeTol()
        )
        {
            convergedD = true;
        }

        if
        (
            (
                solverPerfHE.initialResidual() < solutionTol()
             && residualHE < solutionTol()
            )
         || solverPerfHE.initialResidual() < alternativeTol()
         || residualHE < alternativeTol()
         || absResidualHE < absHETol_
        )
        {
            convergedHE = true;
        }

        if (convergedD && convergedHE)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (he & D), relRes (he & D), matRes, iters (he & D)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfHE.initialResidual()
            << ", " << mag(solverPerfD.initialResidual())
            << ", " << residualHE
            << ", " << residualD
            << ", " << materialResidual
            << ", " << solverPerfHE.nIterations()
            << ", " << solverPerfD.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the enery-momentum loop" << endl;
    }

    return converged;
}


void thermalLinearSolid::readDict()
{
    solidModel::readDict();
    solidModelDict().readIfPresent("absoluteEnergyTolerance", absHETol_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalLinearSolid::thermalLinearSolid(dynamicFvMesh& mesh)
:
    LinearGeomSolid<totalDisplacementSolid>(typeName, mesh),
    absHETol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "absoluteEneryTolerance",
            1e-06
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool thermalLinearSolid::evolve()
{
    Info<< "Evolving thermal solid solver" << endl;

    int iCorr = 0;
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<scalar> solverPerfHE;
    SolverPerformance<vector>::debug = 0;
    SolverPerformance<scalar>::debug = 0;

    Info<< "Solving coupled energy and displacements equation for T and D"
        << endl;

    this->mesh().update();

    blastThermo& thermo = thermal().thermo();
    volScalarField& he = thermo.he();

    // Momentum-energy coupling outer loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        he.storePrevIter();

        // Heat equation
        fvScalarMatrix heEqn
        (
            fvm::ddt(thermo.rho(), he)
          + thermal().divq()
        );

        // Under-relaxation the linear system
        heEqn.relax();

        // Solve the linear system
        solverPerfHE = heEqn.solve();

        // Under-relax the field
        he.relax();

        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(sigma(), "div(sigma)")
          + rho()*g()
          + mechanical().RhieChowCorrection(D(), gradD())
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the field
        relaxField(D(), iCorr);

        // Update increment of displacement
        DD() = D() - D().oldTime();

        // Update velocity
        U() = fvc::ddt(D());

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Update impKf to improve convergence
        // Note: impK and rImpK are not updated as they are used for traction
        // boundaries
        if (iCorr % 10 == 0)
        {
            impKf_ = mechanical().impKf();
        }
    }
    while
    (
        !converged(iCorr, solverPerfD, solverPerfHE, D(), he)
     && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    Info<< "Max T = " << gMax(thermo.T()) << ", "
        << "Min T = " << gMin(thermo.T()) << endl;

    return true;
}


bool thermalLinearSolid::write(const bool write) const
{
    bool good = true;
    if (write)
    {
        const volScalarField& T = thermal().thermo().T();

        Info<< "Max T = " << max(T).value() << nl
            << "Min T = " << min(T).value() << endl;

        // Heat flux
        volVectorField heatFlux
        (
            volVectorField::New
            (
                "heatFlux",
                -this->thermal().k()*fvc::grad(T)
            )
        );
        good = heatFlux.write();

        Info<< "Max magnitude of heat flux = " << max(mag(heatFlux)).value()
            << endl;
    }
    return good && LinearGeomSolid<totalDisplacementSolid>::write(write);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
