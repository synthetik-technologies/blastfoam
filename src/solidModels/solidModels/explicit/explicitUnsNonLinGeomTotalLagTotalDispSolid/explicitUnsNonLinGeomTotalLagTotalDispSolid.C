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

#include "explicitUnsNonLinGeomTotalLagTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "wedgeFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitUnsNonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitUnsNonLinGeomTotalLagTotalDispSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void explicitUnsNonLinGeomTotalLagTotalDispSolid::updateStress()
{
    // Update increment of displacement
    DD() = D() - D().oldTime();

    // Interpolate D to pointD
    mechanical().interpolate(D(), pointD(), false);

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    F_ = I + gradD().T() + tensor::I;
    relF_ = F_ & inv(F_.oldTime());
    J_ = det(F_);
    relJ_ = det(relF_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    impK_ = mechanical().impK();
    impKf_ = mechanical().impKf();
    rImpK_ = 1.0/impK_;

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Increment of displacement
    DD() = D() - D().oldTime();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitUnsNonLinGeomTotalLagTotalDispSolid::explicitUnsNonLinGeomTotalLagTotalDispSolid
(
    dynamicFvMesh& mesh
)
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
    relF_
    (
        IOobject
        (
            "Finv",
            mesh.time().timeName(),
            mesh
        ),
        F_ & inv(F_.oldTime())
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
    relJ_
    (
        IOobject
        (
            "relJ",
            mesh.time().timeName(),
            mesh
        ),
        det(relF_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    // RhieChowScaleFactor_
    // (
    //     solidModelDict().lookupOrDefault<scalar>
    //     (
    //         "RhieChowScale", 0.0
    //     )
    // ),
    JSTScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "JSTScaleFactor", 0.01
        )
    ),
    waveSpeed_
    (
        IOobject
        (
            "waveSpeed",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(Foam::sqrt(impK_/rho()))
    ),
    energies_(mesh, solidModelDict()),
    a_
    (
        IOobject
        (
            "a",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero", dimVelocity/dimTime, vector::zero
        ),
        "zeroGradient"
    )
{
    a_.oldTime();
    U().oldTime();

    // Update stress
    updateStress();

    // Update initial acceleration
    a_.primitiveFieldRef() =
        fvc::div(sigma(), "div(sigma)")().internalField()
       /(rho().internalField());
    a_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitUnsNonLinGeomTotalLagTotalDispSolid::setDeltaT(Time& runTime)
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value
    const scalar requiredDeltaT =
        1.0/
        gMax
        (
            DimensionedField<scalar, Foam::surfaceMesh>
            (
                mesh().surfaceInterpolation::deltaCoeffs().internalField()
               *waveSpeed_.internalField()
            )
        );

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.7071);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    Info<< "maxCo = " << maxCo << nl
        << "deltaT = " << newDeltaT << nl << endl;

    runTime.setDeltaT(newDeltaT);
}


bool explicitUnsNonLinGeomTotalLagTotalDispSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        Info<< "Solving the momentum equation for D" << endl;

        // Central difference scheme

        const dimensionedScalar& deltaT = time().deltaT();
        const dimensionedScalar& deltaT0 = time().deltaT0();

        // Compute the velocity
        // Note: this is the velocity at the middle of the time-step
        U() = U().oldTime() + 0.5*(deltaT + deltaT0)*a_.oldTime();

        // Compute displacement
        D() = D().oldTime() + deltaT*U();

        // Enforce boundary conditions on the displacement field
        D().correctBoundaryConditions();

        // Update the stress field based on the latest D field
        updateStress();

        // Compute acceleration
        // Note the inclusion of a linear bulk viscosity pressure term to
        // dissipate high frequency energies, and a Rhie-Chow term to avoid
        // checker-boarding
        a_.primitiveFieldRef() =
            (
                fvc::div
                (
                    (fvc::interpolate(J_*(sigma() & T(inv(F_)))) & mesh().Sf())
                  + mesh().Sf()*energies_.viscousPressure
                    (
                        rho(), waveSpeed_, gradD()
                    )
                )().primitiveField()
              - JSTScaleFactor_*fvc::laplacian
                (
                    mesh().magSf(),
                    fvc::laplacian
                    (
                        0.5*(deltaT + deltaT0)*impKf_, U(),
                        "laplacian(DU,U)"
                    ),
                    "laplacian(DU,U)"
                )().primitiveField()
            )/rho().primitiveField()
          + g().value();
        a_.correctBoundaryConditions();

        // Check energies
        energies_.checkEnergies
        (
            rho(), U(), D(), DD(), sigma(), gradD(), gradDD(), waveSpeed_, g(),
            0.0, impKf_
        );
    }
    while (mesh().update());

    return true;
}


tmp<vectorField> explicitUnsNonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
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
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField Finv(inv(F_.boundaryField()[patchID]));

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

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
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


scalar explicitUnsNonLinGeomTotalLagTotalDispSolid::CoNum() const
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value
    surfaceScalarField amaxSf(waveSpeed_*mesh().magSf());

    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );
    return 0.5*gMax(sumAmaxSf/mesh().V().field())*mesh().time().deltaTValue();
}


scalar explicitUnsNonLinGeomTotalLagTotalDispSolid::maxCoNum() const
{
    return
        mesh().time().controlDict().lookupOrDefault<scalar>
        (
            "maxCo",
            0.7071
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
