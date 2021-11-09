/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coupledMultiphaseCompressibleSystem.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledMultiphaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        coupledMultiphaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMultiphaseCompressibleSystem::coupledMultiphaseCompressibleSystem
(
    const fvMesh& mesh
)
:
    multiphaseCompressibleSystem(mesh),
    volumeFraction_
    (
        IOobject
        (
            "continuousVolumeFraction",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        0.0
    ),
    alphaRho_
    (
        IOobject
        (
            "alphaRho",
            mesh.time().timeName(),
            mesh
        ),
        volumeFraction_*rho_
    )
{
    forAll(alphas_, phasei)
    {
        volumeFraction_ += alphas_[phasei];
    }
    dynamicCast<multiphaseFluidBlastThermo>
    (
        thermoPtr_()
    ).setTotalVolumeFractionPtr(volumeFraction_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledMultiphaseCompressibleSystem::~coupledMultiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMultiphaseCompressibleSystem::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();

    surfaceScalarField alphaf
    (
        fluxScheme_->interpolate(volumeFraction_, "alpha")
    );

    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) // alphaRhoUPhi_
      - p_*fvc::grad(alphaf)
      - g_*alphaRho_
    );
    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_) // alphaRhoEPhi
      - volumeFraction_*(rhoU_ & g_)
    );

    this->storeAndBlendOld(volumeFraction_);
    volumeFraction_.storePrevIter();

    volumeFraction_ = 0;
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei]) - alphas_[phasei]*fvc::div(phi_)
        );
        this->storeAndBlendDelta(deltaAlpha);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
        this->storeAndBlendDelta(deltaAlphaRho);

        this->storeAndBlendOld(alphas_[phasei]);
        alphas_[phasei] -= dT*deltaAlpha;
        alphas_[phasei].correctBoundaryConditions();
        volumeFraction_ += alphas_[phasei];

        this->storeAndBlendOld(alphaRhos_[phasei]);
        rho_ += alphaRhos_[phasei];
        alphaRhos_[phasei].storePrevIter();
        alphaRhos_[phasei] -= dT*deltaAlphaRho;
    }

    //- Store "old" total density
    rho_ /= max(volumeFraction_.prevIter(), 1e-6);
    rho_.storePrevIter();

    //- Compute new density
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        rho_ += alphaRhos_[phasei];
    }
    rho_ /= max(volumeFraction_, 1e-6);

    thermoPtr_->solve();

    deltaRhoE -= ESource();

    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    this->storeAndBlendOld(rhoU_);
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);

    this->storeAndBlendOld(rhoE_);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::coupledMultiphaseCompressibleSystem::postUpdate()
{
    this->decode();

    compressibleBlastSystem::postUpdate();
}


void Foam::coupledMultiphaseCompressibleSystem::calcAlphas()
{}


void Foam::coupledMultiphaseCompressibleSystem::decode()
{
    volumeFraction_ = min(1.0, max(0.0, 1.0 - alphadPtr_()));
    alphaRho_ = Zero;
    volScalarField sumAlpha
    (
        volScalarField::New
        (
            "sumAlpha",
            rho_.mesh(),
            0.0
        )
    );
    // Correct all but the last phase
    for (label phasei = 0; phasei < alphas_.size() - 1; phasei++)
    {
        alphas_[phasei].maxMin(0.0, 1.0);
        alphas_[phasei].correctBoundaryConditions();
        sumAlpha += alphas_[phasei];

        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());
        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei].boundaryFieldRef() =
            alphas_[phasei].boundaryField()
           *rhos_[phasei].boundaryField();
        alphaRho_ += alphaRhos_[phasei];
    }

    // Correct the last phase
    label lastPhase = alphas_.size() - 1;
    alphas_[lastPhase] = volumeFraction_ - sumAlpha;
    alphas_[lastPhase].maxMin(0.0, 1.0);

    alphaRhos_[lastPhase].max(0.0);
    rhos_[lastPhase] =
            alphaRhos_[lastPhase]
           /max
            (
                alphas_[lastPhase],
                thermo_.thermo(lastPhase).residualAlpha()
            );
    rhos_[lastPhase].correctBoundaryConditions();
    alphaRhos_[lastPhase].boundaryFieldRef() =
        alphas_[lastPhase].boundaryField()
       *rhos_[lastPhase].boundaryField();
    alphaRho_ += alphaRhos_[lastPhase];

    // Update density
    rho_ = alphaRho_/max(volumeFraction_, 1e-10);

    // Update velocity
    volScalarField alphaRhos(alphaRho_);
    alphaRhos.max(1e-10);
    U_.ref() = rhoU_()/alphaRhos();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() =
        alphaRho_.boundaryField()*U_.boundaryField();

    //- Update internal energy
    e_.ref() = rhoE_()/alphaRhos() - 0.5*magSqr(U_());

    thermoPtr_->correct();

    // Update total energy since e may have changed
    rhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}


void Foam::coupledMultiphaseCompressibleSystem::encode()
{
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }
    rho_ = alphaRho_/max(volumeFraction_, 1e-10);
    rhoU_ = alphaRho_*U_;
    rhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}

// ************************************************************************* //
