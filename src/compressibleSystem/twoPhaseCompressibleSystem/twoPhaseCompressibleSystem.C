/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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

#include "twoPhaseCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "MULES.C"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        twoPhaseCompressibleSystem,
        twoPhase
    );
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

void Foam::twoPhaseCompressibleSystem::setModels()
{
    compressibleBlastSystem::setModels();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseCompressibleSystem::twoPhaseCompressibleSystem
(
    const fvMesh& mesh
)
:
    compressibleBlastSystem(mesh, twoPhaseFluidBlastThermo::typeName),
    thermo_
    (
        refCast<twoPhaseFluidBlastThermo>(thermoPtr_())
    ),
    alpha1_(thermo_.alpha1()),
    rho1_(thermo_.thermo(0).rho()),
    alpha2_(thermo_.alpha2()),
    rho2_(thermo_.thermo(1).rho()),
    alphaRho1_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho1_.group()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1_*rho1_
    ),
    alphaRho2_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho2_.group()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha2_*rho2_
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alphaRho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
    ),
    alphaRhoPhi1_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", alphaRho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
    ),
    alphaRhoPhi2_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", alphaRho2_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
    ),
    transportPhaseDensity_(this->lookupOrDefault("transportPhaseDensity", false))
{
    this->fluxScheme_ = fluxScheme::NewMulti(phi_);
    fluxScheme_->phases().insert(alpha1_.group());
    fluxScheme_->phases().insert(alpha2_.group());

    rho_ = alphaRho1_ + alphaRho2_;

    thermo_.initializeModels();
    this->setModels();

    if
    (
        this->lookupOrDefault<bool>("initialDecode", false)
     && alpha1_.headerOk()
     && alphaRho1_.headerOk()
     && alphaRho2_.headerOk()
     && rhoU_.headerOk()
     && rhoE_.headerOk()
    )
    {
        Info<< "Decoding conservative fields"<<endl;
        decode();
    }
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseCompressibleSystem::~twoPhaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseCompressibleSystem::update()
{
    decode();

    // Alpha 1
    autoPtr<ReconstructionScheme<scalar>> alphaLimiter
    (
        ReconstructionScheme<scalar>::New(alpha1_, "alpha", alpha1_.group(), true)
    );
    surfaceScalarField alpha1Own(alphaLimiter->interpolateOwn());
    surfaceScalarField alpha1Nei(alphaLimiter->interpolateNei());

    autoPtr<ReconstructionScheme<scalar>> rho1Limiter
    (
        ReconstructionScheme<scalar>::New(alphaRho1_, "rho", rho1_.group(), true)
    );
    surfaceScalarField alphaRho1Own(rho1Limiter->interpolateOwn());
    surfaceScalarField alphaRho1Nei(rho1Limiter->interpolateNei());

    autoPtr<ReconstructionScheme<scalar>> rho2Limiter
    (
        ReconstructionScheme<scalar>::New(alphaRho2_, "rho", rho2_.group(), true)
    );
    surfaceScalarField alphaRho2Own(rho2Limiter->interpolateOwn());
    surfaceScalarField alphaRho2Nei(rho2Limiter->interpolateNei());


    fluxScheme_->update
    (
        (alphaRho1Own + alphaRho2Own)(),
        (alphaRho1Nei + alphaRho2Nei)(),
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );

    alphaPhi_ = fluxScheme_->flux(alpha1Own, alpha1Nei, phi_);

    // Limit volume fraction flux
    {
        const scalar vDotResidualAlpha = 1e-4;

        // Create copy of alpha to blend and set old time for MULES correction
        volScalarField alpha1Old(alpha1_);
        this->storeAndBlendOld(alpha1Old, false);
        alpha1Old.oldTime();
        alpha1Old.storeOldTimes();

        // Blend fluxes for ODE solver
        this->storeAndBlendDelta(alphaPhi_);

        volScalarField vDot(alpha1_*fvc::div(phi_));
        volScalarField Su
        (
            volScalarField::New
            (
                "Su",
                mesh(),
                dimensionedScalar(vDot.dimensions(), 0.0)
            )
        );
        volScalarField Sp
        (
            volScalarField::New
            (
                "Sp",
                mesh(),
                dimensionedScalar(vDot.dimensions(), 0.0)
            )
        );
        forAll(vDot, celli)
        {
            if (vDot[celli] > 0.0)
            {
                Sp[celli] -=
                    vDot[celli]/max(1.0 - alpha1_[celli], vDotResidualAlpha);
                Su[celli] +=
                    vDot[celli]/max(1.0 - alpha1_[celli], vDotResidualAlpha);
            }
            else if (vDot[celli] < 0.0)
            {
                Sp[celli] += vDot[celli]/max(alpha1_[celli], vDotResidualAlpha);
            }
        }
        this->storeAndBlendDelta(Su);
        this->storeAndBlendDelta(Sp);

        surfaceScalarField phi(phi_);
        this->storeAndBlendDelta(phi);

        // Limit volume fraction flux to ensure boundedness
        MULES::limit
        (
            1.0/mesh().time().deltaT().value(),
            geometricOneField(),
            alpha1Old,
            phi,
            alphaPhi_,
            Sp,
            Su,
            oneField(),
            zeroField(),
            false
        );

        // Using the total field, un-blend volume fraction flux
        alphaPhi_ = this->calcAndStoreDelta(alphaPhi_);
    }

    alphaRhoPhi1_ = fluxScheme_->flux(alphaRho1Own, alphaRho1Nei, phi_);
    alphaRhoPhi2_ = fluxScheme_->flux(alphaRho2Own, alphaRho2Nei, phi_);

    if (mesh().cacheTemporaryObject(alphaRho1Own.name()))
    {
        mesh().cacheTemporaryObject(alphaRho1Own);
        mesh().cacheTemporaryObject(alphaRho1Nei);
    }
    if (mesh().cacheTemporaryObject(alphaRho2Own.name()))
    {
        mesh().cacheTemporaryObject(alphaRho2Own);
        mesh().cacheTemporaryObject(alphaRho2Nei);
    }

    thermo_.update();
}


void Foam::twoPhaseCompressibleSystem::solve()
{
    compressibleBlastSystem::solve();

    // Update changes in volume fraction and phase mass
    volScalarField deltaAlpha
    (
        fvc::div(alphaPhi_) - alpha1_*fvc::div(phi_)
    );
    this->fvTimeInt_->addDeltaSource(alpha1_.name(), deltaAlpha);

    volScalarField deltaAlphaRho1(fvc::div(alphaRhoPhi1_));
    this->fvTimeInt_->addDeltaSource(alphaRho1_.name(), deltaAlphaRho1);

    volScalarField deltaAlphaRho2(fvc::div(alphaRhoPhi2_));
    this->fvTimeInt_->addDeltaSource(alphaRho2_.name(), deltaAlphaRho2);


    // Delta T value
    dimensionedScalar dT = rho_.time().deltaT();

    // Blend old values
    this->storeAndBlendOld(alpha1_, false);
    this->storeAndBlendOld(alphaRho1_);
    alphaRho1_.storePrevIter();

    this->storeAndBlendOld(alphaRho2_);
    alphaRho2_.storePrevIter();

    // Update "old" total density
    rho_ = alphaRho1_ + alphaRho2_;
    rho_.storePrevIter();


    // Blend deltas
    this->storeAndBlendDelta(deltaAlpha);
    this->storeAndBlendDelta(deltaAlphaRho1);
    this->storeAndBlendDelta(deltaAlphaRho2);


    // Update phase 1 volume fraction
    alpha1_ -= dT*deltaAlpha;
    alpha1_.maxMin(0.0, 1.0);
    alpha1_.correctBoundaryConditions();
    alpha2_ = 1.0 - alpha1_;

    // Update phase 1 mass
    alphaRho1_ -= dT*deltaAlphaRho1;
    alphaRho1_.correctBoundaryConditions();

    // Update phase 2 mass
    alphaRho2_ -= dT*deltaAlphaRho2;
    alphaRho2_.correctBoundaryConditions();

    // Update "new" total density
    rho_ = alphaRho1_ + alphaRho2_;


    // Primitive transport of phase densities
    if (transportPhaseDensity_)
    {
        volScalarField divU(fvc::div(phi_));
        volScalarField deltaRho1
        (
            fvc::div(fluxScheme_->flux(rho1_, phi_)) - rho1_*divU
        );
        this->fvTimeInt_->addDeltaSource(rho1_.name(), deltaRho1);

        volScalarField deltaRho2
        (
            fvc::div(fluxScheme_->flux(rho2_, phi_)) - rho2_*divU
        );
        this->fvTimeInt_->addDeltaSource(rho2_.name(), deltaRho2);

        this->storeAndBlendOld(rho1_);
        this->storeAndBlendOld(rho2_);
        this->storeAndBlendDelta(deltaRho1);
        this->storeAndBlendDelta(deltaRho2);

        rho1_ -= dT*deltaRho1;
        rho2_ -= dT*deltaRho2;
    }

    thermo_.solve();
}


void Foam::twoPhaseCompressibleSystem::postUpdate()
{
    this->decode();

    bool updateRho = false;

    // Solve volume fraction
    if (needSolve(alpha1_.name()))
    {
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha1_) - fvc::ddt(alpha1_)
        ==
            models().source(alpha1_)
        );
        constraints().constrain(alphaEqn);
        alphaEqn.solve();
        constraints().constrain(alpha1_);

        alpha1_.maxMin(0.0, 1.0);
        alpha1_.correctBoundaryConditions();
        alpha2_ = 1.0 - alpha1_;

        alphaRho1_ = alpha1_*rho1_;
        alphaRho2_ = alpha2_*rho2_;

        updateRho = true;
    }
    // Solve phase 1 mass
    if (needSolve(rho1_.name()))
    {
        fvScalarMatrix alphaRho1Eqn
        (
            fvm::ddt(alpha1_, rho1_) - fvc::ddt(alphaRho1_)
          + fvm::ddt(thermoPtr_->residualAlpha(), rho1_)
          - fvc::ddt(thermoPtr_->residualAlpha(), rho1_)
        ==
            models().source(alpha1_, rho1_)
        );
        constraints().constrain(alphaRho1Eqn);
        alphaRho1Eqn.solve();
        constraints().constrain(rho1_);

        alphaRho1_ = alpha1_*rho1_;

        updateRho = true;

    }
    // Solve phase 2 mass
    if (needSolve(rho2_.name()))
    {
        fvScalarMatrix alphaRho2Eqn
        (
            fvm::ddt(alpha2_, rho2_) - fvc::ddt(alphaRho2_)
          + fvm::ddt(thermoPtr_->residualAlpha(), rho2_)
          - fvc::ddt(thermoPtr_->residualAlpha(), rho2_)
        ==
            models().source(alpha2_, rho2_)
        );
        constraints().constrain(alphaRho2Eqn);
        alphaRho2Eqn.solve();
        constraints().constrain(rho2_);

        alphaRho2_ = alpha2_*rho2_;

        updateRho = true;
    }

    // Update phase masses
    if (updateRho)
    {
        rho_.storePrevIter();
        rho_ = alphaRho1_ + alphaRho2_;
    }

    compressibleBlastSystem::postUpdate();
}


void Foam::twoPhaseCompressibleSystem::decode()
{
    // Calculate densities
    alphaRho1_.max(0);
    alphaRho2_.max(0);
    const scalar rAlpha1(thermo_.thermo(0).residualAlpha().value());
    const scalar rAlpha2(thermo_.thermo(1).residualAlpha().value());

    if (transportPhaseDensity_)
    {
        // Update densities that have a valid volume fraction
        forAll(alpha1_, celli)
        {
            const scalar alpha1 = alpha1_[celli];
            if (alpha1 > rAlpha1)
            {
                rho1_[celli] = alphaRho1_[celli]/alpha1;
            }
            const scalar alpha2 = alpha2_[celli];
            if (alpha2 > rAlpha2)
            {
                rho2_[celli] = alphaRho2_[celli]/alpha2_[celli];
            }
        }
    }
    else
    {
        rho1_.ref() = alphaRho1_()/max(alpha1_(), rAlpha1);
        rho2_.ref() = alphaRho2_()/max(alpha2_(), rAlpha2);
    }
    rho1_.correctBoundaryConditions();
    rho2_.correctBoundaryConditions();

    alphaRho1_.correctBoundaryConditions();
    alphaRho1_.boundaryFieldRef() = alpha1_.boundaryField()*rho1_.boundaryField();

    alphaRho2_.correctBoundaryConditions();
    alphaRho2_.boundaryFieldRef() = alpha2_.boundaryField()*rho2_.boundaryField();

    rho_ = alphaRho1_ + alphaRho2_;
    compressibleBlastSystem::decode();
}


void Foam::twoPhaseCompressibleSystem::encode()
{
    alphaRho1_ = alpha1_*rho1_;
    alphaRho2_ = alpha2_*rho2_;
    rho_ = alphaRho1_ + alphaRho2_;
    compressibleBlastSystem::encode();
}

// ************************************************************************* //
