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

#include "fluidPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidPhaseModel, 0);
    addToRunTimeSelectionTable(phaseModel, fluidPhaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::fluidPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index,
    const word& thermoType
)
:
    phaseModel(fluid, phaseName, index),
    thermoPtr_
    (
        fluidBlastThermo::New
        (
            fluid.mesh(),
            thermoType,
            this->name_
        )
    ),
    rho_(thermoPtr_->rhoRef()),
    e_(thermoPtr_->he()),
    T_(thermoPtr_->T()),
    p_(thermoPtr_->p()),
    K_
    (
        IOobject
        (
            IOobject::groupName("K", phaseName),
            fluid.mesh().time().name(),
            fluid.mesh()
        ),
        0.5*magSqr(U_)
    ),
    fluxScheme_()
{
    thermoPtr_->read(phaseDict_);
    fluxScheme_ = phaseFluxScheme::New(phi_, residualAlpha().value());

    this->turbulence_ =
        phaseCompressible::momentumTransportModel::New
        (
            *this,
            rho_,
            U_,
            alphaRhoPhi_,
            phi_,
            *this
        );


    fluid.mesh().schemes().setFluxRequired(U_.name());

    this->thermophysicalTransport_ =
        PhaseThermophysicalTransportModel
        <
            phaseCompressible::momentumTransportModel,
            transportThermoModel
        >::New(turbulence_, thermoPtr_());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::~fluidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidPhaseModel::ESource() const
{
    return (*this)*thermoPtr_->ESource();
}


void Foam::fluidPhaseModel::solve()
{
    // Solve momentum and energy transport
    const dimensionedScalar dT = rho().time().deltaT();

    tmp<volVectorField> tgradAlpha = fvc::grad(fluxScheme_->alphaf());
    const volVectorField& gradAlpha = tgradAlpha();
    tmp<volScalarField> deltaAlpha, deltaRho;
    if (solveAlpha_)
    {
        volScalarField& alpha(*this);
        deltaAlpha =
            volScalarField::New
            (
                IOobject::groupName("deltaAlpha", name_),
                // (fluid_.U() & gradAlpha)
                fvc::div(fluid_.phi()*fluxScheme_->alphaf())
              - alpha*fvc::div(fluid_.phi())
              + fluxScheme_->alphaCorrector(*this)
            );
        this->fvTimeInt_->addDeltaSource(alpha.name(), deltaAlpha.ref());

        // deltaRho =
        //     volScalarField::New
        //     (
        //         IOobject::groupName("deltaRho", name_),
        //         fvc::div(fluxScheme_->flux(rho_, phi_))
        //       - rho_*fvc::div(phi_)
        //     );
    }

    volScalarField deltaAlphaRho
    (
        IOobject::groupName("deltaAlphaRho", name_),
        fvc::div(alphaRhoPhi_)
    );
    this->fvTimeInt_->addDeltaSource(alphaRho_.name(), deltaAlphaRho);

    volVectorField deltaAlphaRhoU
    (
        IOobject::groupName("deltaAlphaRhoU", name_),
        fvc::div(alphaRhoUPhi_)
      - fluid_.PI()*gradAlpha
      - alphaRho_*fluid_.g()
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

    if (fluid_.hasMassTransfer(*this))
    {
        forAll(fluid_.phases(), phasei)
        {
            const phaseModel& otherPhase = fluid_.phases()[phasei];
            if (&otherPhase != this && fluid_.hasMassTransfer(*this, otherPhase))
            {
                volScalarField mD(fluid_.mDot(*this, otherPhase));
                volScalarField alphaD(fluid_.mDotByRho(mD, *this, otherPhase));
                if (solveAlpha_)
                {
                    deltaAlpha.ref() -= alphaD;
                }
                deltaAlphaRho -= mD;
                deltaAlphaRhoU -= fluid_.mDotU(mD, *this, otherPhase);
                deltaAlphaRhoE -=
                    fluid_.mDotE(mD, *this, otherPhase);
            }
        }
    }

    if (explicitViscosity_ && turbulence_.valid())
    {
        tmp<volSymmTensorField> tdevTau(turbulence_->devTau());
        // tdevTau.ref() += (2.0/3.0)*rhoEff()*turbulence_->k()*symmTensor::I;

        deltaAlphaRhoU += fvc::div(tdevTau());
        deltaAlphaRhoE +=
            fvc::div
            (
                fvc::dotInterpolate(mesh().Sf(), tdevTau)
              & flux().Uf()
            )
          + fvc::div(thermophysicalTransport_->q()*mesh().magSf());
    }


    // Transport volume fraction if required
    if (solveAlpha_)
    {
        this->storeAndBlendDelta(deltaAlpha.ref());

        volScalarField& alpha(*this);
        this->storeAndBlendOld(alpha);
        alpha.storePrevIter();

        alpha -= dT*deltaAlpha;
        alpha.max(0.0);
        // alpha.maxMin(0.0, 1.0);
        alpha.correctBoundaryConditions();

        // // Estimate rho
        // this->storeAndBlendDelta(deltaRho.ref());
        // this->storeAndBlendOld(rho_);
        // rho_ -= dT*deltaRho;
    }

    this->storeAndBlendOld(alphaRho_);
    this->storeAndBlendDelta(deltaAlphaRho);
    alphaRho_.storePrevIter();
    alphaRho_ -= this->mesh().time().deltaT()*deltaAlphaRho;
    alphaRho_.max(0);

    this->storeAndBlendOld(alphaRhoU_);
    alphaRhoU_.storePrevIter();
    this->storeAndBlendDelta(deltaAlphaRhoU);
    alphaRhoU_ -= cmptMultiply(dT*deltaAlphaRhoU, solutionDs_);

    this->storeAndBlendOld(alphaRhoE_);
    alphaRhoE_.storePrevIter();
    this->storeAndBlendDelta(deltaAlphaRhoE);
    alphaRhoE_ -= dT*deltaAlphaRhoE;

    // Solve thermodynamic models (activation and afterburn)
    thermoPtr_->solve();
}


void Foam::fluidPhaseModel::solveExplicit()
{}


void Foam::fluidPhaseModel::storeExplicit()
{
    thermoPtr_->storeExplicit();

    if
    (
        solveAlpha_
     && needSolve(static_cast<const volScalarField&>(*this).name())
    )
    {
        alphaAdvection_ = fvc::ddt(*this);
    }

    if (needSolve(rho().name()))
    {
        alphaRhoAdvection_ = fvc::ddt(alphaRho_);
    }

    if (needSolve(U_.name()) || turbulence_.valid())
    {
        alphaRhoUAdvection_ = fvc::ddt(alphaRhoU_);
    }

    if (needSolve(he().name()) || turbulence_.valid())
    {
        alphaRhoEAdvection_ = fvc::ddt(alphaRhoE_);
    }
}


void Foam::fluidPhaseModel::solveImplicit()
{
    volScalarField& alpha(*this);
    if (needSolve(alpha.name()) && solveAlpha_)
    {
        //- Solve momentum equation (implicit stresses)
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha) - alphaAdvection_()
         ==
            models().source(alpha)
        );
        constraints().constrain(alphaEqn);
        alphaEqn.solve();
        constraints().constrain(alpha);
    }

    if (needSolve(rho().name()))
    {
        //- Solve momentum equation (implicit stresses)
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(alpha, rho()) - alphaRhoAdvection_()
          + fvm::ddt(this->residualAlpha(), rho())
          - fvc::ddt(this->residualAlpha(), rho())
         ==
            models().source(alpha, rho())
        );
        constraints().constrain(rhoEqn);
        rhoEqn.solve();
        constraints().constrain(rho());

        alphaRho_ = alpha*rho();
    }

    // Viscous
    if (turbulence_.valid())
    {
        turbulence_->predict();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->predict();
    }

    tmp<surfaceVectorField> devTau;
    if (needSolve(U_.name()) || (turbulence_.valid() && !explicitViscosity_))
    {
        tmp<fvVectorMatrix> divDevTau;
        if (!explicitViscosity_ && turbulence_.valid())
        {
            divDevTau = turbulence_->divDevTau(U_);
        }

        fvVectorMatrix UEqn
        (
            fvm::ddt(alpha, rho(), U_) - alphaRhoUAdvection_()
          + fvc::ddt(this->residualAlphaRho(), U_)
          - fvm::ddt(this->residualAlphaRho(), U_)
         ==
            models().source(*this, rho(), U_)
        );

        if (divDevTau.valid())
        {
            UEqn += divDevTau();
        }

        UEqn.relax();
        constraints().constrain(UEqn);
        UEqn.boundaryManipulate(U_.boundaryFieldRef());

        UEqn.solve();
        constraints().constrain(U_);

        if (divDevTau.valid())
        {
            devTau = divDevTau().flux();
        }

        K_ = 0.5*magSqr(U_);
        alphaRhoU_ = alphaRho_*U_;
    }

    // Solve thermal energy diffusion
    if (needSolve(he().name()) || turbulence_.valid())
    {
        fvScalarMatrix EEqn
        (
            fvm::ddt(alpha, rho(), he()) - alphaRhoEAdvection_()
          + fvc::ddt(alpha, rho(), K_)
          + fvc::ddt(this->residualAlphaRho(), he())
          - fvm::ddt(this->residualAlphaRho(), he())
         ==
            models().source(*this, rho(), he())
        );

        if (devTau.valid())
        {
            EEqn +=
                fvc::div(devTau & flux().Uf())
              + thermophysicalTransport_->divq(he());
        }

        EEqn.relax();

        constraints().constrain(EEqn);
        EEqn.boundaryManipulate(he().boundaryFieldRef());
        EEqn.solve();
        constraints().constrain(he());

        alphaRhoE_ = alphaRho_*(he() + K_);
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->correct();
    }

    thermo().solveImplicit();
    thermo().correct();
}


void Foam::fluidPhaseModel::postUpdate()
{
    thermoPtr_->postUpdate();
}

void Foam::fluidPhaseModel::clear()
{
    fluxScheme_->clear();
    alphaAdvection_.clear();
    alphaRhoAdvection_.clear();
    alphaRhoUAdvection_.clear();
    alphaRhoEAdvection_.clear();
    thermoPtr_->clear();
}


void Foam::fluidPhaseModel::update()
{
    const volScalarField& alpha = *this;
    {
        fluxScheme_->update
        (
            alpha,
            rho_,
            U_,
            e_,
            p_,
            thermoPtr_->speedOfSound(),
            phi_,
            alphaRhoPhi_,
            alphaRhoUPhi_,
            alphaRhoEPhi_
        );
    }

    thermoPtr_->update();
    phaseModel::update();
}


void Foam::fluidPhaseModel::correctDeltas()
{
    const dimensionedScalar& deltaT = this->time().deltaT();
    if (solveAlpha_)
    {
        volScalarField& alpha(*this);
        tmp<volScalarField> deltaAlpha
        (
            volScalarField::New
            (
                IOobject::groupName("deltaAlpha", name_),
                (alpha - alpha.prevIter())/deltaT
            )
        );
        calcAndStoreDelta(deltaAlpha());
    }

    tmp<volVectorField> deltaAlphaRhoU
    (
        volVectorField::New
        (
            IOobject::groupName("deltaAlphaRhoU", name_),
            (alphaRhoU_ - alphaRhoU_.prevIter())/deltaT
        )
    );
    calcAndStoreDelta(deltaAlphaRhoU());

    tmp<volScalarField> deltaAlphaRhoE
    (
        volScalarField::New
        (
            IOobject::groupName("deltaAlphaRhoE", name_),
            (alphaRhoE_ - alphaRhoE_.prevIter())/deltaT
        )
    );
    calcAndStoreDelta(deltaAlphaRhoE());
}

void Foam::fluidPhaseModel::decode()
{
    const fvConstraints& constraints = this->constraints();

    this->correctBoundaryConditions();
    const volScalarField& alpha = *this;
    volScalarField alphaLimited(Foam::max(*this, residualAlpha()));
    const scalar rAlpha(residualAlpha().value());

    // Update density (only if alpha > residual)
    forAll(*this, celli)
    {
        const scalar alpha = (*this)[celli];
        if (alpha > rAlpha)
        {
            rho_[celli] = alphaRho_[celli]/alpha;
        }
    }
    // rho_.internalFieldRef() = alphaRho_()/alphaLimited();
    phaseFluxScheme::correctPhaseFields(alpha, rho_, rAlpha);
    rho_.max(thermo().residualRho());
    rho_.correctBoundaryConditions();

    alphaRho_ == (*this)*rho_;
    // alphaRho_.correctBoundaryConditions();
    // alphaRho_.boundaryFieldRef() ==
    //     (*this).boundaryField()*rho_.boundaryField();
    volScalarField alphaRhoLimited(Foam::max(alphaRho_, residualAlphaRho()));


    // Update velocity
    U_.internalFieldRef() = alphaRhoU_()/(alphaRhoLimited());
    phaseFluxScheme::correctPhaseFields(alpha, U_, rAlpha);
    if (constraints.constrainsField(U_.name()))
    {
        constraints.constrain(U_);
    }
    constraints.constrain(U_);
    U_.correctBoundaryConditions();

    K_ = 0.5*magSqr(U_);

    alphaRhoU_ == alphaRho_*U_;
    // alphaRhoU_.correctBoundaryConditions();
    // alphaRhoU_.boundaryFieldRef() ==
    //     (*this).boundaryField()*rho_.boundaryField()*U_.boundaryField();


    // Update internal energy
    const scalar rAlphaRho = this->residualAlphaRho().value();
    forAll(e_, celli)
    {
        const scalar alphaRho = alphaRho_[celli];
        if (alphaRho > rAlphaRho)
        {
            e_[celli] = alphaRhoE_[celli]/alphaRho - K_[celli];
        }
        // else
        // {
        //     e_[celli] = thermoPtr_->cellhe(thermoPtr_->TLow(), celli);
        // }
    }
    // e_.internalFieldRef() = alphaRhoE_()/alphaRhoLimited() - K_();
    phaseFluxScheme::correctPhaseFields(alpha, e_, rAlpha);
    constraints.constrain(e_);
    e_.correctBoundaryConditions();

    // Correct thermo (T, mu, alpha, etc)
    thermoPtr_->correct();

    // Update total energy because e may have changed
    alphaRhoE_ == alphaRho_*(e_ + K_);

    // Limit speed of sound in cells with low volume fraction
    volScalarField& c = thermoPtr_->speedOfSound();
    forAll(c, celli)
    {
        if (alpha[celli] < rAlpha || alphaRho_[celli] < rAlphaRho)
        {
            c[celli] = small;
            p_[celli] = small;
        }
    }
    phaseFluxScheme::correctPhaseFields(alpha, c, rAlpha);
    c.correctBoundaryConditions();

    // phaseFluxScheme::correctPhaseFields(alpha, p_, rAlpha);
    constraints.constrain(p_);
    p_.correctBoundaryConditions();
}


void Foam::fluidPhaseModel::encode()
{
    K_ = 0.5*magSqr(U_);

    alphaRho_ = (*this)*rho_;
    alphaRhoU_ = alphaRho_*U_;
    alphaRhoE_ = alphaRho_*(e_ + K_);
}


Foam::tmp<Foam::volVectorField>
Foam::fluidPhaseModel::gradP() const
{
    return fvc::grad(fluxScheme_->pf());
}


Foam::tmp<Foam::volVectorField>
Foam::fluidPhaseModel::gradAlpha() const
{
    return fvc::grad(fluxScheme_->alphaf());
}



void Foam::fluidPhaseModel::correctThermo()
{
    thermoPtr_->correct();
}

// ************************************************************************* //
