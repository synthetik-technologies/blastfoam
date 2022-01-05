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

#include "nonLinGeomTotalLagSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagSolid::nonLinGeomTotalLagSolid(dynamicFvMesh& mesh)
:
    solidModel(typeName, mesh, nonLinGeom(), incremental()),
    F_
    (
        IOobject
        (
            "F",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{
    DDisRequired();

    //- Dummy Call to make sure the necessary old fields are initialized
    fvc::d2dt2(rho().oldTime(), D().oldTime());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    SolverPerformance<vector> solverPerfDD;
    SolverPerformance<vector>::debug = 0;


    Info<< "Solving the total Lagrangian form of the momentum equation for DD"
        << endl;

    // Reset enforceLinear switch
    enforceLinear() = false;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Momentum equation incremental displacement total Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho(), DD())
          + fvc::d2dt2(rho().oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
          + rho()*g()
          + stabilisation().stabilisation(DD(), gradDD(), impK_)
        );

        // Enforce linear to improve convergence
        if (enforceLinear())
        {
            // Replace nonlinear terms with linear
            // Note: the mechanical law could still be nonlinear
            DDEqn +=
                fvc::div(J_*Finv_ & sigma(), "div(sigma)")
              - fvc::div(sigma());
        }

        // Under-relax the linear system
        DDEqn.relax();

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field using fixed or adaptive under-relaxation
        relaxField(DD(), iCorr);

        // Update the total displacement
        D() = D().oldTime() + DD();

        // Update gradient of displacement increment
        mechanical().grad(DD(), gradDD());

        // Update the gradient of total displacement
        gradD() = gradD().oldTime() + gradDD();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Check if outer loops are diverging
        if (!enforceLinear())
        {
            checkEnforceLinear(J_);
        }

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DDEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Update impKf to improve convergence
        // Note: impK and rImpK are not updated as they are used for traction
        // boundaries
        // if (iCorr % 10 == 0)
        // {
        //     impKf_ = mechanical().impKf();
        // }
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

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD(), pointDD());

    // Increment of displacement
    D() = D().oldTime() + DD();

    // Increment of point displacement
    pointD() = pointD().oldTime() + pointDD();

    // Velocity
    U() = fvc::ddt(D());

    return true;
}


tmp<vectorField> nonLinGeomTotalLagSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    if (enforceLinear())
    {
        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - n*pressure)
                  - (n & pSigma)
                  + impK*(n & pGradDD)
                )*rImpK
            )
        );
    }
    else
    {
        // Patch total deformation gradient inverse
        const tensorField& Finv = Finv_.boundaryField()[patchID];

        // Patch unit normals (deformed configuration)
        vectorField nCurrent(Finv.T() & n);
        nCurrent /= mag(nCurrent);

        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - nCurrent*pressure)
                  - (nCurrent & pSigma)
                  + impK*(n & pGradDD)
                )*rImpK
            )
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
