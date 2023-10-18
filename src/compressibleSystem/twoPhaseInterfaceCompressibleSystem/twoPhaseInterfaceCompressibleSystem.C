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

#include "twoPhaseInterfaceCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "EulerDdtScheme.H"
#include "gaussConvectionScheme.H"
#include "MULES.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseInterfaceCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        twoPhaseInterfaceCompressibleSystem,
        twoPhase
    );
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::twoPhaseInterfaceCompressibleSystem::rhoUSource() const
{
    return
        compressibleBlastSystem::rhoUSource()
      + fvc::reconstruct(interface_.surfaceTensionForce()*mesh().magSf());
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseInterfaceCompressibleSystem::rhoESource() const
{
    return
        compressibleBlastSystem::rhoESource()
      + (
            fvc::reconstruct(interface_.surfaceTensionForce()*mesh().magSf())
          & U()
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseInterfaceCompressibleSystem::twoPhaseInterfaceCompressibleSystem
(
    const fvMesh& mesh
)
:
    twoPhaseCompressibleSystem(mesh),
    interfaceSystem(U_, *this),
    interface_(alpha1_, alpha2_, U_, *this),
    psi_
    (
        IOobject
        (
            IOobject::groupName("levelSet", alpha1_.group()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->levelSet(alpha1_)
    ),
    nHatf_
    (
        IOobject
        (
            "nHatf",
            mesh.time().timeName(),
            mesh
        ),
        (fvc::interpolate(fvc::grad(psi_)) & mesh_.Sf())
       /(mag(fvc::interpolate(fvc::grad(psi_))) + small)
    )
    // surfaceTension_(surfaceTensionModel::New(*this, mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseInterfaceCompressibleSystem::~twoPhaseInterfaceCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseInterfaceCompressibleSystem::update()
{
    // Decode primitives
    decode();

    // tmp<surfaceScalarField> talpha1Own, talpha1Nei;
    // autoPtr<ReconstructionScheme<scalar>> alpha1Limiter
    // (
    //     ReconstructionScheme<scalar>::New(alpha1_, "alpha", alpha1_.group(), true)
    // );
    tmp<surfaceScalarField> talpha1Own = fvc::interpolate
    (
        alpha1_,
        phi_,
        reconstruction::scheme("alpha", alpha1_.group(), mesh(), true, true)
    );
    tmp<surfaceScalarField> talpha1Nei(talpha1Own());
    // alpha1Limiter->interpolateOwnNei(talpha1Own, talpha1Nei);
    const surfaceScalarField& alpha1Own = talpha1Own();
    const surfaceScalarField& alpha1Nei = talpha1Nei();

    surfaceScalarField alpha2Own(1.0 - alpha1Own);
    surfaceScalarField alpha2Nei(1.0 - alpha1Nei);


    tmp<surfaceScalarField> trho1Own, trho1Nei;
    autoPtr<ReconstructionScheme<scalar>> rho1Limiter
    (
        ReconstructionScheme<scalar>::New(rho1_, "rho", rho1_.group(), true)
    );
    rho1Limiter->interpolateOwnNei(trho1Own, trho1Nei);

    tmp<surfaceScalarField> trho2Own, trho2Nei;
    autoPtr<ReconstructionScheme<scalar>> rho2Limiter
    (
        ReconstructionScheme<scalar>::New(rho2_, "rho", rho2_.group(), true)
    );
    rho2Limiter->interpolateOwnNei(trho2Own, trho2Nei);

    // Adjust densities on phase boundaries so that densities with zero volume
    // fraction are not used
    if (!transportPhaseDensity_)
    {
        fluxScheme::correctPhaseFields
        (
            alpha1_,
            trho1Own.ref(), trho1Nei.ref(),
            thermo_.thermo(0).residualAlpha().value()
        );
        fluxScheme::correctPhaseFields
        (
            alpha2_,
            trho2Own.ref(), trho2Nei.ref(),
            thermo_.thermo(1).residualAlpha().value()
        );
    }

    // Compute total phase masses
    tmp<surfaceScalarField> talphaRho1Own = surfaceScalarField::New
    (
        rho1Limiter->ownName(alphaRho1_.name()),
        alpha1Own*trho1Own()
    );
    tmp<surfaceScalarField> talphaRho1Nei = surfaceScalarField::New
    (
        rho1Limiter->neiName(alphaRho1_.name()),
        alpha1Nei*trho1Nei()
    );
    if (mesh().cacheTemporaryObject(talphaRho1Own().name()))
    {
        mesh().cacheTemporaryObject(talphaRho1Own.ref());
        mesh().cacheTemporaryObject(talphaRho1Nei.ref());
    }

    tmp<surfaceScalarField> talphaRho2Own = surfaceScalarField::New
    (
        rho2Limiter->ownName(alphaRho2_.name()),
        alpha2Own*trho2Own()
    );
    tmp<surfaceScalarField> talphaRho2Nei = surfaceScalarField::New
    (
        rho2Limiter->ownName(alphaRho2_.name()),
        alpha2Nei*trho2Nei()
    );
    if (mesh().cacheTemporaryObject(talphaRho2Own().name()))
    {
        mesh().cacheTemporaryObject(talphaRho2Own.ref());
        mesh().cacheTemporaryObject(talphaRho2Nei.ref());
    }

    surfaceScalarField rhoOwn
    (
        reconstruction::ownName(rho_.name()),
        talphaRho1Own() + talphaRho2Own()
    );
    surfaceScalarField rhoNei
    (
        reconstruction::neiName(rho_.name()),
        talphaRho1Nei() + talphaRho2Nei()
    );
    if (mesh().cacheTemporaryObject(rhoOwn.name()))
    {
        mesh().cacheTemporaryObject(rhoOwn);
        mesh().cacheTemporaryObject(rhoNei);
    }

    // Compute fluxes using the reconstructed total density
    fluxScheme_->update
    (
        rhoOwn,
        rhoNei,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );

    // Interface fields
    interface_.correct();

    // Limit volume fraction flux
    {
        // Create copy of alpha to blend and set old time for MULES correction
        volScalarField alpha1Old(alpha1_);
        this->storeAndBlendOld(alpha1Old, false);
        alpha1Old.oldTime();
        alpha1Old.storeOldTimes();

        // Blend fluxes for ODE solver
        alphaPhi_ = fluxScheme_->flux(alpha1Own, alpha1Nei, phi_);
        this->storeAndBlendDelta(alphaPhi_);

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
            zeroField(),
            zeroField(),//(alpha1_.v()*fvc::div(phi_)().v())(),
            oneField(),
            zeroField(),
            false
        );

        // Using the total field, un-blend volume fraction flux
        alphaPhi_ = this->calcAndStoreDelta(alphaPhi_);
    }

    // Update phase mass fluxes
    alphaRhoPhi1_ = fluxScheme_->phaseFlux
    (
        alpha1_,
        trho1Own(), trho1Nei(),
        alphaPhi_,
        thermo_.thermo(0).residualAlpha().value()
    );
    surfaceScalarField alphaPhi2(phi_ - alphaPhi_);
    alphaRhoPhi2_ = fluxScheme_->phaseFlux
    (
        alpha2_,
        trho2Own(), trho2Nei(),
        alphaPhi2,
        thermo_.thermo(1).residualAlpha().value()
    );

    // UPtrList<volScalarField> alphas(2);
    // alphas.set(0, &alpha1_);
    // alphas.set(1, &alpha2_);
    //
    // UPtrList<volScalarField> rhos(2);
    // rhos.set(0, &rho1_);
    // rhos.set(1, &rho2_);
    //
    // UPtrList<surfaceScalarField> alphaRhoPhis(2);
    // alphaRhoPhis.set(0, &alphaRhoPhi1_);
    // alphaRhoPhis.set(1, &alphaRhoPhi2_);
    //
    // limitAlphaRhoPhis(alphas, rhos, alphaRhoPhis, phi_, rhoPhi_);

    // Update thermo
    thermo_.update();
}


// void Foam::twoPhaseInterfaceCompressibleSystem::decode()
// {
//     // Limit first phase volume fraction and calculate second phase volume fraction
//     alpha1_.maxMin(0.0, 1.0);
//     alpha1_.correctBoundaryConditions();
//     alpha2_ = 1.0 - alpha1_;
//
//     // Limit phase masses
//     alphaRho1_.max(0);
//     alphaRho2_.max(0);
//
//
//     // Only update cells that have a valid volume fraction
//     // other cell densities are handled by transport of density
//     const scalar rAlpha1 = thermo_.thermo(0).residualAlpha().value();
//     const scalar rAlpha2 = thermo_.thermo(1).residualAlpha().value();
//     // forAll(alpha1_, celli)
//     // {
//     //     const scalar alpha1 = alpha1_[celli];
//     //     if (alpha1 > rAlpha1)
//     //     {
//     //         rho1_[celli] = alphaRho1_[celli]/alpha1;
//     //     }
//     //     const scalar alpha2 = alpha2_[celli];
//     //     if (alpha2 > rAlpha2)
//     //     {
//     //         rho2_[celli] = alphaRho2_[celli]/alpha2_[celli];
//     //     }
//     // }
//
//     rho1_.ref() = alphaRho1_()/max(alpha1_(), rAlpha1);
//     rho2_.ref() = alphaRho2_()/max(alpha2_(), rAlpha2);
//
//     // const extendedNLevelCFCCellToCellStencil& stencil
//     // (
//     //     extendedNLevelCFCCellToCellStencil::New(mesh(), 3)
//     // );
//     // List<List<scalar>> alpha1Nei(alpha1_.size());
//     // List<List<scalar>> rho1Nei(rho1_.size());
//     // List<List<scalar>> rho2Nei(rho2_.size());
//     //
//     // stencil.collectData(alpha1_, alpha1Nei);
//     // stencil.collectData(rho1_, rho1Nei);
//     // stencil.collectData(rho2_, rho2Nei);
//
//     // forAll(alpha1_, celli)
//     // {
//     //     const scalar alpha1 = alpha1_[celli];
//     //     if (alpha1 < 0.5)
//     //     {
//     //         scalar sumRhoW = 0.0;
//     //         scalar sumW = 0.0;
//     //         forAll(rho1Nei[celli], cj)
//     //         {
//     //             if (alpha1Nei[celli][cj] > 0.5)
//     //             {
//     //                 scalar w = 1.0/max(1.0 - alpha1Nei[celli][cj], rAlpha1);
//     //                 sumRhoW += rho1Nei[celli][cj]*w;
//     //                 sumW += w;
//     //             }
//     //         }
//     //         rho1_[celli] = sumRhoW/max(sumW, rAlpha1);
//     //         // alpha1_[celli] = alphaRho1_[celli]/max(rho1_[celli], rRho1);
//     //     }
//     //
//     //     const scalar alpha2 = alpha2_[celli];
//     //     if (alpha2 < 0.5)
//     //     {
//     //         scalar sumRhoW = 0.0;
//     //         scalar sumW = 0.0;
//     //         forAll(rho2Nei[celli], cj)
//     //         {
//     //             if (alpha1Nei[celli][cj] < 0.5)
//     //             {
//     //                 scalar w = 1.0/max(alpha1Nei[celli][cj], rAlpha2);
//     //                 sumRhoW += rho2Nei[celli][cj]*w;
//     //                 sumW += w;
//     //             }
//     //         }
//     //         rho2_[celli] = sumRhoW/max(sumW, rAlpha1);
//     //         // alpha2_[celli] = alphaRho2_[celli]/max(rho2_[celli], rRho2);
//     //     }
//     // }
//
//
//     // {
//     //     // Gradient of level set function
//     //     surfaceVectorField gradAlphaf(fvc::interpolate(fvc::grad(alpha1_)));
//     //
//     //     // Face unit interface normal
//     //     nHatf_ = (gradAlphaf & mesh_.Sf())/(mag(gradAlphaf) + deltaN_);
//     // }
//     // {
//     //     psi_ = this->levelSet(alpha1_);
//     //
//     //     // Gradient of level set function
//     //     surfaceVectorField gradPsif(fvc::interpolate(fvc::grad(psi_)));
//     //
//     //     // Face unit interface normal
//     //     nHatf_ = (gradPsif & mesh_.Sf())/max(mag(gradPsif), 1e-6);
//     // }
//
//     // this->GFMExtension(rho1_, -psi_);
//     // this->GFMExtension(rho2_, psi_, nHatf_, fvc::div(nHatf_)());
//
//     rho1_.correctBoundaryConditions();
//     rho2_.correctBoundaryConditions();
//     alphaRho1_.boundaryFieldRef() = alpha1_.boundaryField()*rho1_.boundaryField();
//     alphaRho2_.boundaryFieldRef() = alpha2_.boundaryField()*rho2_.boundaryField();
//
//     rho_ = alphaRho1_ + alphaRho2_;
//     compressibleBlastSystem::decode();
// }


// ************************************************************************* //
