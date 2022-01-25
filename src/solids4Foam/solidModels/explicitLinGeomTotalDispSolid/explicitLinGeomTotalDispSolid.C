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

#include "explicitLinGeomTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "wedgeFvPatch.H"

#include "operations.H"
#include "mechanics.H"
#include "gradientSchemes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitLinGeomTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitLinGeomTotalDispSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void explicitLinGeomTotalDispSolid::updateStress()
{
    // Update increment of displacement
    DD() = D() - D().oldTime();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitLinGeomTotalDispSolid::explicitLinGeomTotalDispSolid
(
    dynamicFvMesh& mesh
)
:
    solidModel(typeName, mesh, nonLinGeom(), incremental()),
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
            IOobject::READ_IF_PRESENT,
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
    DisRequired();

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

void explicitLinGeomTotalDispSolid::setDeltaT(Time& runTime)
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


bool explicitLinGeomTotalDispSolid::evolve()
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

        volTensorField P(sigma() & tensor::I);
        volVectorField rhoU(rho()*U());

        operations ops(mesh());
        volVectorField Px(ops.decomposeTensorX(P));
        volVectorField Py(ops.decomposeTensorY(P));
        volVectorField Pz(ops.decomposeTensorZ(P));

        volTensorField gradRhoU(fvc::grad(rhoU));
        volTensorField gradPx(fvc::grad(Px));
        volTensorField gradPy(fvc::grad(Py));
        volTensorField gradPz(fvc::grad(Pz));

        // Reconstruction
        surfaceVectorField rhoUOwn
        (
            surfaceVectorField::New
            (
                "rhoUOwn",
                mesh(),
                dimDensity*dimVelocity
            )
        );
        surfaceVectorField rhoUNei
        (
            surfaceVectorField::New
            (
                "rhoUNei",
                mesh(),
                dimDensity*dimVelocity
            )
        );

        surfaceTensorField POwn
        (
            surfaceTensorField::New
            (
                "POwn",
                mesh(),
                P.dimensions()
            )
        );
        surfaceTensorField PNei
        (
            surfaceTensorField::New
            (
                "PNei",
                mesh(),
                P.dimensions()
            )
        );

        gradientSchemes gradSchemes(D());
        gradSchemes.reconstruct(rhoU, gradRhoU,rhoUOwn, rhoUNei);
        gradSchemes.reconstruct(P, gradPx, gradPy, gradPz, POwn, PNei);

        surfaceVectorField N(mesh().Sf()/mesh().magSf());
        surfaceVectorField tractionOwn(POwn & N);
        surfaceVectorField tractionNei(PNei & N);

        volScalarField pWaveSpeed
        (
            sqrt(mechanical().elasticModulus()/rho())
        );
        volScalarField sWaveSpeed
        (
            sqrt(mechanical().shearModulus()/rho())
        );
        volTensorField F
        (
            volTensorField::New("F", mesh(), tensor::I)
        );
        mechanics mech(F, ops);
        mech.correct(pWaveSpeed, sWaveSpeed);

        // Acoustic Riemann solver
        surfaceVectorField tractionC
        (
            0.5*(tractionOwn + tractionNei)
          + 0.5*(mech.stabRhoU() & (rhoUNei - rhoUOwn))
        );
//         rhoUC =
//             0.5*(rhoUOwn_ + rhoUNei_)
//           + 0.5*(mech.stabTraction() & (tractionNei_ - tractionOwn_));


        // Compute acceleration
        // Note the inclusion of a linear bulk viscosity pressure term to
        // dissipate high frequency energies, and a Rhie-Chow term to avoid
        // checker-boarding
        a_.primitiveFieldRef() =
            (
//                 fvc::surfaceIntegrate(mesh().magSf()*tractionC)()() +
                fvc::div
                (
                    (mesh().Sf() & fvc::interpolate(sigma())) +
                    mesh().Sf()*energies_.viscousPressure
                    (
                        rho(), waveSpeed_, gradD()
                    )
                )().primitiveField()
                // This corresponds to Lax–Friedrichs smoothing
//                 + JSTScaleFactor_*fvc::laplacian
//                   (
//                       0.5*(deltaT + deltaT0)*impKf_,
//                       U(),
//                       "laplacian(DU,U)"
//                   )().primitiveField()
              - JSTScaleFactor_*fvc::laplacian
                (
                    mesh().magSf(),
                    fvc::laplacian
                    (
                        0.5*(deltaT + deltaT0)*impKf_, U(), "laplacian(DU,U)"
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


tmp<vectorField> explicitLinGeomTotalDispSolid::tractionBoundarySnGrad
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


scalar explicitLinGeomTotalDispSolid::CoNum() const
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


scalar explicitLinGeomTotalDispSolid::maxCoNum() const
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
