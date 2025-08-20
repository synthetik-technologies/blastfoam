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
    fluxScheme_(phaseFluxScheme::New(phi_))
{
    thermoPtr_->read(phaseDict_);

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

    tmp<volVectorField> tdelta = fvc::grad(fluxScheme_->deltaAlphaf());
    const volVectorField& delta = tdelta();
    tmp<volScalarField> deltaAlpha;
    if (solveAlpha_)
    {
        volScalarField& alpha(*this);
        deltaAlpha =
            volScalarField::New
            (
                IOobject::groupName("deltaAlpha", name_),
                (fluid_.U() & delta)
              + fluxScheme_->alphaCorrector(*this)
            );
        this->fvTimeInt_->addDeltaSource(alpha.name(), deltaAlpha.ref());
    }

    volVectorField deltaAlphaRhoU
    (
        IOobject::groupName("deltaAlphaRhoU", name_),
        fvc::div(alphaRhoUPhi_)
      - fluid_.PI()*delta
      - alphaRho_*fluid_.g()
    );
    this->fvTimeInt_->addDeltaSource(alphaRhoU_.name(), deltaAlphaRhoU);

    volScalarField deltaAlphaRhoE
    (
        IOobject::groupName("deltaAlphaRhoE", name_),
        fvc::div(alphaRhoEPhi_)
      - ESource()
      - fluid_.PI()*(fluid_.U() & delta)
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
                volScalarField alphaD(fluid_.mDotByRho(*this, otherPhase));
                if (solveAlpha_)
                {
                    deltaAlpha.ref() -= alphaD;
                }
                deltaAlphaRhoU -= fluid_.mDotU(mD, *this, otherPhase);
                deltaAlphaRhoE -=
                    fluid_.mDotE(mD, *this, otherPhase)
                  - alphaD*p();
            }
        }
    }
    this->storeAndBlendDelta(deltaAlphaRhoU);
    this->storeAndBlendDelta(deltaAlphaRhoE);


    this->storeAndBlendOld(alphaRhoU_);
    alphaRhoU_ -= cmptMultiply(dT*deltaAlphaRhoU, solutionDs_);
    alphaRhoU_.correctBoundaryConditions();

    this->storeAndBlendOld(alphaRhoE_);
    alphaRhoE_ -= dT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();

    // Solve phase density transport to store old time values
    phaseModel::solveAlphaRho();

    // Transport volume fraction if required
    if (solveAlpha_)
    {
        this->storeAndBlendDelta(deltaAlpha.ref());

        volScalarField& alpha(*this);
        this->storeAndBlendOld(alpha);

        volScalarField alphaOld(alpha);
        alpha -= dT*deltaAlpha;
        alpha.maxMin(0.0, 1.0);
        alpha.correctBoundaryConditions();
    }

    // Solve thermodynamic models (activation and afterburn)
    thermoPtr_->solve();
}


void Foam::fluidPhaseModel::postUpdate()
{
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
    if (needSolve(rho().name()))
    {
        //- Solve momentum equation (implicit stresses)
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(alpha, rho()) - fvc::ddt(alphaRho_)
          + fvm::ddt(residualAlpha(), rho())
          - fvc::ddt(residualAlpha(), rho())
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


void Foam::fluidPhaseModel::update()
{
    const volScalarField& alpha = *this;
    {
        // Info<<"no solve alpha "<<this->name()<<endl;
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
            alphaRhoEPhi_,
            residualAlpha().value()
        );
    }

    thermoPtr_->update();
    phaseModel::update();
}


void Foam::fluidPhaseModel::decode()
{
    const fvConstraints& constraints = this->constraints();

    this->correctBoundaryConditions();
    const volScalarField& alpha = *this;
    volScalarField alphaLimited(Foam::max(*this, residualAlpha()));

    const scalar rAlpha(residualAlpha().value());
    forAll(*this, celli)
    {
        const scalar alpha = (*this)[celli];
        if (alpha > rAlpha)
        {
            rho_[celli] = alphaRho_[celli]/alpha;
        }
    }
    // rho_.internalFieldRef() = alphaRho_()/alphaLimited();
    rho_.max(thermo().residualRho());
    rho_.correctBoundaryConditions();

    alphaRho_.correctBoundaryConditions();
    alphaRho_.boundaryFieldRef() ==
        (*this).boundaryField()*rho_.boundaryField();
    volScalarField alphaRhoLimited(Foam::max(alphaRho_, residualAlphaRho()));

    U_.internalFieldRef() = alphaRhoU_()/(alphaRhoLimited());
    if (constraints.constrainsField(U_.name()))
    {
        constraints.constrain(U_);
    }
    constraints.constrain(U_);
    U_.correctBoundaryConditions();

    alphaRhoU_.internalFieldRef() = (*this)()*rho_()*U_;
    alphaRhoU_.boundaryFieldRef() ==
        (*this).boundaryField()*rho_.boundaryField()*U_.boundaryField();

    // forAll(e_, celli)
    // {
    //     if ((*this)[celli] > residualAlpha().value())
    //     {
    //         e_[celli] =
    //             alphaRhoE_[celli]/alphaRhoLimited[celli]
    //           - 0.5*magSqr(U_[celli]);
    //     }
    //     else
    //     {
    //         e_[celli] = small;
    //     }
    // }
    e_.internalFieldRef() = alphaRhoE_()/alphaRhoLimited() - 0.5*magSqr(U_());
    constraints.constrain(e_);
    e_.correctBoundaryConditions();

    thermoPtr_->correct();
    thermoPtr_->speedOfSound() *= pos(alpha - residualAlpha());
    thermoPtr_->speedOfSound().max(small);

    // Update total energy because e may have changed
    alphaRhoE_ == alphaRho_*(e_ + 0.5*magSqr(U_));

    constraints.constrain(p_);
}


void Foam::fluidPhaseModel::encode()
{
    alphaRho_ = (*this)*rho_;
    alphaRhoU_ = alphaRho_*U_;
    alphaRhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
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
