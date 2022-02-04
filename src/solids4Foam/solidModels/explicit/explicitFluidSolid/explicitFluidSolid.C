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

#include "explicitFluidSolid.H"
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

defineTypeNameAndDebug(explicitFluidSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitFluidSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void explicitFluidSolid::decode()
{
    U_.ref() = rhoU_()/rho_();

    DD().ref() = U_()*solidModel::mesh().time().deltaT();
    DD().correctBoundaryConditions();

    // Update the total displacement
    D() = D().oldTime() + DD();

    // Interpolate DD to pointDD
    mechanical().interpolate(DD(), pointDD(), false);

    // Update gradient of displacement increment
    mechanical().grad(DD(), pointDD(), gradDD());

    // Update the gradient of total displacement
    gradD() = gradD().oldTime() + gradDD();

    // Relative deformation gradient
    relF_ = I + gradDD().T();

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();


    U_.boundaryFieldRef() =
        DD().boundaryField()/solidModel::mesh().time().deltaTValue();

    rhoU_.boundaryFieldRef() =
        rho_.boundaryField()*U_.boundaryField();

    e_.ref() = rhoE_()/rho_() - 0.5*magSqr(U_());

    thermo_.correct();

    //- Update total energy because the e field may have been modified
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));

    mechanical().correct(sigma());
}


void explicitFluidSolid::encode()
{
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitFluidSolid::explicitFluidSolid
(
    dynamicFvMesh& mesh
)
:
    updatedLagSolid<incrementalSolid>(typeName, mesh, false),
    timeIntegrationSystem(typeName, mesh),
    integrator_(timeIntegrator::NewRef(mesh)),
    thermo_
    (
        dynamicCast<fluidBlastThermo>(thermal().thermo())
    ),
    rho_(thermo_.rho()),
    p_(thermo_.p()),
    T_(thermo_.T()),
    e_(thermo_.he()),
    U_(U()),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimDensity*dimVelocity, Zero),
        "zeroGradient"
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity*sqr(dimVelocity), 0.0)
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity*dimVelocity*dimArea, 0.0)
    ),
    rhoUPhi_
    (
        IOobject
        (
            "rhoUPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("0", dimDensity*sqr(dimVelocity)*dimArea, Zero)
    ),
    rhoEPhi_
    (
        IOobject
        (
            "rhoEPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity*pow3(dimVelocity)*dimArea, 0.0)
    ),
    fluxScheme_(fluxScheme::NewSingle(mesh)),
    LFScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "LFScaleFactor", 0.001
        )
    ),
    JSTScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "JSTScaleFactor", 0.01
        )
    ),
    energies_(mesh, solidModelDict())
{
    integrator_.addSystem(*this);
    encode();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitFluidSolid::update()
{
    decode();
    fluxScheme_->update
    (
        rho_,
        U_,
        thermo_.he(),
        thermo_.p(),
        thermo_.speedOfSound(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
    thermo_.update();
}


void explicitFluidSolid::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();

    volScalarField deltaRho("deltaRho", fvc::div(rhoPhi_));

    surfaceVectorField sigmaFlux
    (
        (fluxScheme_->interpolate(sigma()) & solidModel::mesh().Sf())
      + energies_.viscousPressure
        (
            rho(),
            fvc::interpolate(thermo_.speedOfSound()),
            gradD()
        )*solidModel::mesh().Sf()
    );

    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_)
      - fvc::div(sigmaFlux)
      - g()*rho_
    );

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_)
      - fvc::div(sigmaFlux & fluxScheme_->Uf())
      - thermo_.ESource()
      - (rhoU_ & g())
    );

//     if (LFScaleFactor_ > small)
//     {
//         volVectorField LFSmoothing
//         (
//             LFScaleFactor_*fvc::laplacian
//             (
//                 dT*impKf_,
//                 U_,
//                 "laplacian(DU,U)"
//             )
//         );
//         deltaRhoU += LFSmoothing;
//         deltaRhoE += LFSmoothing & U_;
//     }
//     if (JSTScaleFactor_ > small)
//     {
//         volVectorField JSTSmoothing
//         (
//             JSTScaleFactor_*fvc::laplacian
//             (
//                 solidModel::mesh().magSf(),
//                 fvc::laplacian
//                 (
//                     dT*impKf_,
//                     U_,
//                     "laplacian(DU,U)"
//                 ),
//                 "laplacian(DU,U)"
//             )
//         );
//         deltaRhoU -= JSTSmoothing;
//         deltaRhoE -= JSTSmoothing & U_;
//     }

    //- Store old values
    this->storeAndBlendOld(rho_);
    rho_.storePrevIter();
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRho);
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    //- Solve for momentum and energy
    rho_ -= dT*deltaRho;
    rho_.correctBoundaryConditions();

    rhoU_ -= dT*deltaRhoU;
    rhoE_ -= dT*deltaRhoE;

    thermo_.solve();
}


void explicitFluidSolid::postUpdate()
{}


void explicitFluidSolid::setDeltaT(Time& runTime)
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
        1.0
        /max
        (
            solidModel::mesh().surfaceInterpolation::deltaCoeffs()
           *fvc::interpolate(thermo_.speedOfSound())
        ).value();

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.7071);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    Info<< "maxCo = " << maxCo << nl
        << "deltaT = " << newDeltaT << nl << endl;

    runTime.setDeltaT(newDeltaT);
}


bool explicitFluidSolid::evolve()
{

    integrator_.preUpdateMesh();

    // Mesh update loop
//     do
    {
        Info<< "Solving the momentum equation for D" << endl;

        integrator_.integrate();

        surfaceScalarField waveSpeed
        (
            fvc::interpolate(thermo_.speedOfSound())
        );

        // Check energies
        energies_.checkEnergies
        (
            rho_,
            U_,
            D(),
            DD(),
            sigma(),
            gradD(),
            gradDD(),
            waveSpeed,
            g(),
            0.0,
            impKf_
        );

        integrator_.clear();
    }
//     while (solidModel::mesh().update());

    return true;
}


scalar explicitFluidSolid::CoNum() const
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value

    surfaceScalarField amaxSf
    (
        fvc::interpolate(thermo_.speedOfSound())*solidModel::mesh().magSf()
    );

    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(solidModel::mesh().boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    amaxSf += mag(phi_);

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );
    return
        0.5*gMax
        (
            sumAmaxSf/solidModel::mesh().V().field()
        )*solidModel::mesh().time().deltaTValue();
}


scalar explicitFluidSolid::maxCoNum() const
{
    return
        solidModel::mesh().time().controlDict().lookupOrDefault<scalar>
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
