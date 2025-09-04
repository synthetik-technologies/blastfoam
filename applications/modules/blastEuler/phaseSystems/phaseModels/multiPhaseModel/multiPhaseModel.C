/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "multiPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "MULES.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiPhaseModel, 0);
    addToRunTimeSelectionTable
    (
        phaseModel,
        multiPhaseModel,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * * //

void Foam::multiPhaseModel::updateFluxes
(
    const PtrList<surfaceScalarField>& alphasOwn,
    const PtrList<surfaceScalarField>& alphasNei,
    const PtrList<surfaceScalarField>& alphaRhosOwn,
    const PtrList<surfaceScalarField>& alphaRhosNei
)
{
    surfaceScalarField alphaOwn
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                reconstruction::ownName("alpha"),
                name_
            ),
            mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );
    surfaceScalarField alphaNei
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                reconstruction::neiName("alpha"),
                name_
            ),
            mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );
    surfaceScalarField rhoOwn
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                reconstruction::ownName("rho"),
                name_
            ),
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    surfaceScalarField rhoNei
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                reconstruction::neiName("rho"),
                name_
            ),
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    forAll(alphaRhoPhis_, phasei)
    {
        alphaOwn += alphasOwn[phasei];
        alphaNei += alphasNei[phasei];

        rhoOwn += alphaRhosOwn[phasei];
        rhoNei += alphaRhosNei[phasei];
    }

    fluxScheme_->update
    (
        *this,
        alphaOwn,
        alphaNei,
        rho_,
        rhoOwn,
        rhoNei,
        U_,
        e_,
        p_,
        speedOfSound(),
        phi_,
        alphaRhoPhi_,
        alphaRhoUPhi_,
        alphaRhoEPhi_
    );

    // fluxScheme_->update
    // (
    //     rhoOwn,
    //     rhoNei,
    //     U_,
    //     e_,
    //     p_,
    //     speedOfSound()(),
    //     phi_,
    //     rhoPhi_,
    //     rhoUPhi_,
    //     rhoEPhi_
    // );

    forAll(alphaRhoPhis_, phasei)
    {
        alphaPhis_[phasei] = fluxScheme_->flux
        (
            alphasOwn[phasei],
            alphasNei[phasei],
            phi_
        );
        alphaRhoPhis_[phasei] = fluxScheme_->flux
        (
            alphaRhosOwn[phasei],
            alphaRhosNei[phasei],
            phi_
        );
    }
/*
    if (MUSLESLimiting_)
    {
        tmp<volScalarField> tdivPhi(fvc::div(phi_));
        const volScalarField& divPhi = tdivPhi();
        // Limit alpha flux
        surfaceScalarField phi(phi_);
        this->storeAndBlendDelta(phi);

        UPtrList<const volScalarField> alphas(alphas_.size());
        forAll(alphas_, phasei)
        {
            const volScalarField& alpha = alphas_[phasei];
            volScalarField Su
            (
                IOobject::groupName("Su", alpha.group()),
                alphas_[phasei]*divPhi
            );
            this->storeAndBlendDelta(Su);

            surfaceScalarField& alphaPhi = alphaPhis_[phasei];
            this->storeAndBlendDelta(alphaPhi);

            alphas.set(phasei, &alphas_[phasei]);
            volScalarField alphaOld(alphas_[phasei]);
            this->blendOld(alphaOld);
            alphaOld.storeOldTimes();

            MULES::limit
            (
                1.0/mesh().time().deltaT().value(),
                geometricOneField(),
                alphaOld,
                phi,
                alphaPhi,
                zeroField(),
                Su,
                oneField(),
                zeroField(),
                false
            );
            alphaPhi = this->calcAndStoreDelta(alphaPhi);
        }
        MULES::limitSum(alphas, alphaPhis_, phi_);
    }*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiPhaseModel::multiPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    fluidPhaseModel
    (
        fluid,
        phaseName,
        index,
        multiphaseFluidBlastThermo::typeName
    ),
    thermo_(dynamicCast<multiphaseFluidBlastThermo>(thermoPtr_())),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size()),
    densityReconstruction_
    (
        phaseDict_.lookupOrDefault("densityReconstruction", false)
    )
{
    thermo_.setTotalVolumeFractionPtr(*this);

    //- Temporarily Store read density
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            fluid.mesh().time().name(),
            fluid.mesh()
        ),
        fluid.mesh(),
        0.0
    );
    alphaRho_ = dimensionedScalar(dimDensity, 0.0);

    wordList phaseNames(alphas_.size());
    forAll(alphas_, phasei)
    {
        phaseNames[phasei] = IOobject::groupName
        (
            thermo_.phaseNames()[phasei],
            this->name()
        );
        sumAlpha += alphas_[phasei];
        word phaseName = phaseNames[phasei];
        alphaRhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", phaseName),
                    fluid.mesh().time().name(),
                    fluid.mesh()
                ),
                alphas_[phasei]*rhos_[phasei],
                rhos_[phasei].boundaryField().types()
            )
        );
        alphaRho_ += alphaRhos_[phasei];
        alphaPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaPhi", phaseName),
                    fluid.mesh().time().name(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
            )
        );
        alphaRhoPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoPhi", phaseName),
                    fluid.mesh().time().name(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
            )
        );
    }

    this->fluxScheme_->phases() = phaseNames;

    // Reset density to correct value
    volScalarField& alpha = *this;
    alpha = sumAlpha;
    alpha.correctBoundaryConditions();

    rho_ = alphaRho_/Foam::max(sumAlpha, residualAlpha());

    solveAlpha(true);

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiPhaseModel::~multiPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiPhaseModel::solve()
{
    if (fluid_.hasMassTransfer(*this))
    {
        NotImplemented;
    }

    // Solve momentum and energy transport
    const dimensionedScalar dT = rho().time().deltaT();

    tmp<volVectorField> tgradAlpha = this->gradAlpha();
    const volVectorField& gradAlpha = tgradAlpha();


    volVectorField deltaAlphaRhoU
    (
        IOobject::groupName("deltaAlphaRhoU", name_),
        fvc::div(alphaRhoUPhi_)
      - fluid_.PI()*gradAlpha
      - (*this)*rho()*fluid_.g() // alphaRho has already been updated
    );
    this->fvTimeInt_->addDeltaSource(alphaRhoU_.name(), deltaAlphaRhoU);

    volScalarField deltaAlphaRhoE
    (
        IOobject::groupName("deltaAlphaRhoE", name_),
        fvc::div(alphaRhoEPhi_)
      - ESource()
      - fluid_.PI()*(fluid_.U() & gradAlpha)
      - (alphaRhoU_ & fluid_.g())
    );
    this->fvTimeInt_->addDeltaSource(alphaRhoE_.name(), deltaAlphaRhoE);

    // if (fluid_.hasMassTransfer(*this))
    // {
    //     forAll(fluid_.phases(), phasei)
    //     {
    //         const phaseModel& otherPhase = fluid_.phases()[phasei];
    //         if (&otherPhase != this && fluid_.hasMassTransfer(*this, otherPhase))
    //         {
    //             volScalarField mD(fluid_.mDot(*this, otherPhase));
    //             volScalarField alphaD(fluid_.mDotByRho(*this, otherPhase));
    //             deltaAlphaRhoU -= fluid_.mDotU(mD, *this, otherPhase);
    //             deltaAlphaRhoE -=
    //                 fluid_.mDotE(mD, *this, otherPhase)
    //               - alphaD*p();
    //         }
    //     }
    // }
    this->storeAndBlendDelta(deltaAlphaRhoU);
    this->storeAndBlendDelta(deltaAlphaRhoE);


    this->storeAndBlendOld(alphaRhoU_);
    alphaRhoU_ -= cmptMultiply(dT*deltaAlphaRhoU, solutionDs_);
    alphaRhoU_.correctBoundaryConditions();

    this->storeAndBlendOld(alphaRhoE_);
    alphaRhoE_ -= dT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();


    dynamicCast<volScalarField>(*this) = 0.0;
    forAll(alphas_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei])
          - alphas_[phasei]*fvc::div(phi_)
          + fluxScheme_->alphaCorrector(alphas_[phasei])
        );
        this->fvTimeInt_->addDeltaSource(alphas_[phasei].name(), deltaAlpha);
        this->storeAndBlendDelta(deltaAlpha);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
        this->fvTimeInt_->addDeltaSource
        (
            alphaRhos_[phasei].name(),
            deltaAlphaRho
        );
        this->storeAndBlendDelta(deltaAlphaRho);

        this->storeAndBlendOld(alphas_[phasei]);
        alphas_[phasei] -= dT*deltaAlpha;
        alphas_[phasei].correctBoundaryConditions();

        this->storeAndBlendOld(alphaRhos_[phasei]);
        alphaRhos_[phasei].storePrevIter();
        alphaRhos_[phasei] -= dT*deltaAlphaRho;
        alphaRhos_[phasei].correctBoundaryConditions();

        *this += alphas_[phasei];
    }


    // Solve thermodynamic models (activation and afterburn)
    thermoPtr_->solve();
}


void Foam::multiPhaseModel::postUpdate()
{
    // Viscous
    if (turbulence_.valid())
    {
        turbulence_->predict();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->predict();
    }

    // Solve phase mass
    bool needUpdate = false;
    alphaRho_.storePrevIter();
    forAll(rhos_, phasei)
    {
        bool alphaRhoUpdate = false;
        volScalarField& alpha(alphas_[phasei]);
        if (needSolve(alpha.name()))
        {
            //- Solve momentum equation (implicit stresses)
            fvScalarMatrix alphaEqn
            (
                fvm::ddt(alpha) - fvc::ddt(alpha)
             ==
                models().source(alpha)
            );
            constraints().constrain(alphaEqn);
            alphaEqn.solve();
            constraints().constrain(alpha);

            alphaRhoUpdate = true;
        }

        volScalarField& rho(rhos_[phasei]);
        if (needSolve(rho.name()))
        {
            dimensionedScalar rAlpha
            (
                thermo_.thermo(phasei).residualAlpha()
            );
            //- Solve momentum equation (implicit stresses)
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(alpha, rho) - fvc::ddt(alphaRhos_[phasei])
              + fvm::ddt(rAlpha, rho)
              - fvc::ddt(rAlpha, rho)
             ==
                models().source(alpha, rho)
            );
            constraints().constrain(rhoEqn);
            rhoEqn.solve();
            constraints().constrain(rho);

            alphaRhoUpdate = true;
        }
        if (alphaRhoUpdate)
        {
            alphaRhos_[phasei] = alpha*rho;
            needUpdate = true;
        }
    }
    if (needUpdate)
    {
        alphaRho_ = alphaRhos_[0];
        for (label phasei = 1; phasei < alphaRhos_.size(); phasei++)
        {
            alphaRho_ += alphaRhos_[phasei];
        }
        rho_ = alphaRho_/Foam::max(*this, residualAlpha());
    }


    volScalarField& alpha(*this);
    if (needSolve(alpha.name()) && solveAlpha_)
    {
        //- Solve momentum equation (implicit stresses)
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha) - fvc::ddt(alpha)
         ==
            models().source(alpha)
        );
        constraints().constrain(alphaEqn);
        alphaEqn.solve();
        constraints().constrain(alpha);
    }

    alphaRho_.storePrevIter();

    dimensionedScalar smallAlphaRho(residualAlphaRho());
    if (needSolve(U_.name()) || turbulence_.valid())
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(alphaRho_, U_) - fvc::ddt(alphaRhoU_)
          + fvc::ddt(smallAlphaRho, U_) - fvm::ddt(smallAlphaRho, U_)
         ==
            models().source(*this, rho(), U_)
        );
        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
            alphaRhoE_ +=
                rho().time().deltaT()
               *fvc::div
                (
                    fvc::dotInterpolate
                    (
                        rho().mesh().Sf(),
                        turbulence_->devTau()
                    )
                  & flux().Uf()
                );
        }
        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        alphaRhoU_ = alphaRho_*U_;

        he() = alphaRhoE_/Foam::max(alphaRho_, smallAlphaRho) - 0.5*magSqr(U_);
    }

    // Solve thermal energy diffusion
    if (needSolve(he().name()) || turbulence_.valid())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(alphaRho_, he())
          - fvc::ddt(alphaRho_.prevIter(), he())
          + fvc::ddt(smallAlphaRho, he())
          - fvm::ddt(smallAlphaRho, he())
         ==
            models().source(*this, rho(), he())
        );

        if (turbulence_.valid())
        {
            // Add thermal energy diffusion
            eEqn += thermophysicalTransport_->divq(he());
        }
        constraints().constrain(eEqn);
        eEqn.solve();
        constraints().constrain(he());

        alphaRhoE_ = alphaRho_*(he() + 0.5*magSqr(U_));
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->correct();
    }

    thermo().postUpdate();
    thermo().correct();
}


void Foam::multiPhaseModel::update()
{
    PtrList<surfaceScalarField> alphasOwn(alphas_.size());
    PtrList<surfaceScalarField> alphasNei(alphas_.size());
    PtrList<surfaceScalarField> alphaRhosOwn(alphas_.size());
    PtrList<surfaceScalarField> alphaRhosNei(alphas_.size());

    phi_ = fvc::relative(fvc::flux(U_), U_);

    forAll(alphas_, phasei)
    {
        autoPtr<ReconstructionScheme<scalar>> alphaLimiter
        (
            ReconstructionScheme<scalar>::New
            (
                alphas_[phasei],
                "alpha",
                alphas_[phasei].group(),
                true
            )
        );
        if (alphaLimiter->upwind())
        {
            alphasOwn.set
            (
                phasei,
                surfaceScalarField::New
                (
                    alphaLimiter->ownName(alphas_[phasei].name()),
                    alphaLimiter->interpolate(phi_)
                )
            );
            alphasNei.set
            (
                phasei,
                surfaceScalarField::New
                (
                    alphaLimiter->neiName(alphas_[phasei].name()),
                    alphasOwn[phasei]
                )
            );
        }
        else
        {
            alphasOwn.set(phasei, alphaLimiter->interpolateOwn());
            alphasNei.set(phasei, alphaLimiter->interpolateNei());
        }
    }

    // Ensure sum of alphaOwn and alphaNei = 1
    {
        forAll(alphasOwn[0], facei)
        {
            scalar sumAlphaOwn = 0;
            scalar sumAlphaNei = 0;
            forAll(alphasOwn, phasei)
            {
                sumAlphaOwn += alphasOwn[phasei][facei];
                sumAlphaNei += alphasNei[phasei][facei];
            }
            forAll(alphasOwn, phasei)
            {
                alphasOwn[phasei][facei] /= sumAlphaOwn;
                alphasNei[phasei][facei] /= sumAlphaNei;
            }
        }

        UPtrList<surfaceScalarField::Boundary> balphasOwn(alphas_.size());
        UPtrList<surfaceScalarField::Boundary> balphasNei(alphas_.size());
        forAll(alphasOwn, phasei)
        {
            balphasOwn.set(phasei, &alphasOwn[phasei].boundaryFieldRef());
            balphasNei.set(phasei, &alphasNei[phasei].boundaryFieldRef());
        }

        forAll(balphasOwn[0], patchi)
        {
            forAll(balphasOwn[0][patchi], facei)
            {
                scalar sumAlphaOwn = 0;
                scalar sumAlphaNei = 0;
                forAll(balphasOwn, phasei)
                {
                    sumAlphaOwn += balphasOwn[phasei][patchi][facei];
                    sumAlphaNei += balphasNei[phasei][patchi][facei];
                }
                forAll(alphasOwn, phasei)
                {
                    balphasOwn[phasei][patchi][facei] /= sumAlphaOwn;
                    balphasNei[phasei][patchi][facei] /= sumAlphaNei;
                }
            }
        }
    }

    forAll(alphas_, phasei)
    {
        if (densityReconstruction_)
        {
            autoPtr<ReconstructionScheme<scalar>> rhoLimiter
            (
                ReconstructionScheme<scalar>::New
                (
                    rhos_[phasei],
                    "rho",
                    rhos_[phasei].group(),
                    true
                )
            );
            surfaceScalarField rhoOwn(rhoLimiter->interpolateOwn());
            surfaceScalarField rhoNei(rhoLimiter->interpolateNei());
            // if (!transportPhaseDensity_)
            {
                fluxScheme::correctPhaseFields
                (
                    alphas_[phasei],
                    rhos_[phasei],
                    rhoOwn, rhoNei,
                    thermo_.thermo(phasei).residualAlpha().value()
                );
            }
            alphaRhosOwn.set
            (
                phasei,
                surfaceScalarField::New
                (
                    rhoLimiter->ownName(alphaRhos_[phasei].name()),
                    alphasOwn[phasei]*rhoOwn
                )
            );
            alphaRhosNei.set
            (
                phasei,
                surfaceScalarField::New
                (
                    rhoLimiter->neiName(alphaRhos_[phasei].name()),
                    alphasNei[phasei]*rhoNei
                )
            );
        }
        else
        {
            autoPtr<ReconstructionScheme<scalar>> rhoLimiter
            (
                ReconstructionScheme<scalar>::New
                (
                    alphaRhos_[phasei],
                    "rho",
                    rhos_[phasei].group(),
                    true
                )
            );
            alphaRhosOwn.set(phasei, rhoLimiter->interpolateOwn());
            alphaRhosNei.set(phasei, rhoLimiter->interpolateNei());
        }
    }

    updateFluxes(alphasOwn, alphasNei, alphaRhosOwn, alphaRhosNei);
    phaseModel::update();
    thermoPtr_->update();
}


void Foam::multiPhaseModel::scaleVolumeFraction
(
    const scalar sumAlpha,
    const label celli
)
{
    forAll(alphas_, phasei)
    {
        alphas_[phasei][celli] /= sumAlpha;
    }
    (*this)[celli] /= sumAlpha;
}


void Foam::multiPhaseModel::correctVolumeFraction
(
    const scalar alpha,
    const label celli
)
{
    scalar sumAlpha = ::Foam::max((*this)[celli], residualAlpha().value());
    forAll(alphas_, phasei)
    {
        alphas_[phasei][celli] *= alpha/sumAlpha;
    }


    (*this)[celli] = alpha;
}


void Foam::multiPhaseModel::decode()
{
    const fvConstraints& constraints = this->constraints();

    // Calculate densities
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);

    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /Foam::max
            (
                alphas_[phasei],
                thermo_.thermo(phasei).residualAlpha()
            );
        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];

        alphaRho_ += alphaRhos_[phasei];
    }
    volScalarField& alpha = *this;
    this->correctBoundaryConditions();

    rho_ = alphaRho_/Foam::max(alpha, residualAlpha());

    volScalarField alphaRhoLimited(alphaRho_);
    alphaRhoLimited.max(1e-10);
    U_.internalFieldRef() = alphaRhoU_()/(alphaRhoLimited());
    if (constraints.constrainsField(U_.name()))
    {
        constraints.constrain(U_);
        alphaRhoU_.internalFieldRef() = (*this)()*rho_()*U_;
    }
    U_.correctBoundaryConditions();

    alphaRhoU_.correctBoundaryConditions();
    alphaRhoU_.boundaryFieldRef() ==
        alphaRho_.boundaryField()*U_.boundaryField();

    e_.internalFieldRef() = alphaRhoE_()/alphaRhoLimited() - 0.5*magSqr(U_());
    constraints.constrain(e_);
    e_.correctBoundaryConditions();

    thermoPtr_->correct();
    thermoPtr_->speedOfSound() *= pos(alpha - residualAlpha());
    thermoPtr_->speedOfSound().max(small);

    // Update total energy because e may have changed
    alphaRhoE_ == alphaRho_*(e_ + 0.5*magSqr(U_));
}


void Foam::multiPhaseModel::encode()
{
    //- Scale volume fractions to new value
    alphaRho_ = dimensionedScalar(dimDensity, 0.0);
    volScalarField& alpha(*this);
    alpha = 0.0;
    forAll(alphas_, phasei)
    {
        alpha += alphas_[phasei];
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }

    alphaRhoU_ == alphaRho_*U_;
    alphaRhoE_ == alphaRho_*(e_ + 0.5*magSqr(U_));
}


// ************************************************************************* //
