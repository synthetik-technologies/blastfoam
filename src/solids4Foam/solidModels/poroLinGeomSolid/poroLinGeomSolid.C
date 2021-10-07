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

#include "poroLinGeomSolid.H"
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

defineTypeNameAndDebug(poroLinGeomSolid, 0);
addToRunTimeSelectionTable(solidModel, poroLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool poroLinGeomSolid::converged
(
    const int iCorr,
    const SolverPerformance<vector>& solverPerfD,
    const SolverPerformance<scalar>& solverPerfp,
    const volVectorField& D,
    const volScalarField& p
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate displacement residual
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

    // Calculate pressure residual
    const scalar residualp =
        gMax
        (
            DimensionedField<double, volMesh>
            (
                mag(p.internalField() - p.prevIter().internalField())
               /max
                (
                    gMax
                    (
                        DimensionedField<double, volMesh>
                        (
                            mag(p.internalField() - p.oldTime().internalField())
                        )
                    ),
                    SMALL
                )
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        if
        (
            mag(solverPerfD.initialResidual()) < solutionTol()
         && solverPerfp.initialResidual() < solutionTol()
         && residualD < solutionTol()
         && residualp < solutionTol()
        )
        {
            Info<< "    All residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            residualD < alternativeTol()
         && residualp < alternativeTol()
        )
        {
            Info<< "    The relative residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            mag(solverPerfD.initialResidual()) < alternativeTol()
         && solverPerfp.initialResidual() < alternativeTol()
        )
        {
            Info<< "    The solver residuals have converged" << endl;
            converged = true;
        }
        else
        {
            converged = false;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, resD, resP, relResD, relResP, matRes, iters" << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << mag(solverPerfD.initialResidual())
            << ", " << solverPerfp.initialResidual()
            << ", " << residualD
            << ", " << residualp
            << ", " << materialResidual
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
            << "Max iterations reached within momentum-pressure loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

poroLinGeomSolid::poroLinGeomSolid(dynamicFvMesh& mesh)
:
    solidModel(typeName, mesh, nonLinGeom(), incremental()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
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
    gradp_
    (
        IOobject
        (
            "grad(p)",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimPressure/dimLength, vector::zero)
    ),
    hydraulicConductivity_(solidModelDict().lookup("hydraulicConductivity")),
    gammaWater_(solidModelDict().lookup("waterSpecificWeight")),
    porosity_(solidModelDict().lookup("porosity")),
    saturation_(solidModelDict().lookup("degreeOfSaturation")),
    KWater_(solidModelDict().lookup("waterBulkModulus")),
    rKprime_
    (
        (saturation_/KWater_)
      + (1.0 - saturation_)
        /dimensionedScalar("atmosphericPressure", dimPressure, 1e+05)
    )
{
    // Store old time of p
    p_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool poroLinGeomSolid::evolve()
{
    Info << "Evolving poro solid solver" << endl;

    int iCorr = 0;
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<scalar> solverPerfp;
    SolverPerformance<vector>::debug = 0;

    Info<< "Solving the pressure equation for p and momentum equation for D"
        << endl;

    // Pressure-displacement coupling outer loop
    do
    {
        // Pressure equation

        // Store fields for under-relaxation and residual calculation
        p_.storePrevIter();

        // Pressure equation
        fvScalarMatrix pEqn
        (
            (porosity_*rKprime_)*fvm::ddt(p_)
          + fvc::div(U())
         == (hydraulicConductivity_/gammaWater_)*fvm::laplacian(p_)
        );

        // Under-relaxation the linear system
        pEqn.relax();

        // Solve the linear system
        solverPerfp = pEqn.solve();

        // Under-relax the field
        p_.relax();

        // Update gradient of pressure
        gradp_ = fvc::grad(p_);


        // Momentum equation

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

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the field
        relaxField(D(), iCorr);

        // Update increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Update velocity as it is used in the pEqn
        U() = fvc::ddt(D());

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
        !converged(iCorr, solverPerfD, solverPerfp, D(), p_)
     && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    return true;
}


tmp<vectorField> poroLinGeomSolid::tractionBoundarySnGrad
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
