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

#include "granularPhaseModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(granularPhaseModel, 0);
    addToRunTimeSelectionTable(phaseModel, granularPhaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularPhaseModel::granularPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    phaseModel(fluid, phaseName, index),
    kineticTheoryModel
    (
        *this,
        phaseDict_.subDict("kineticTheoryCoeffs")
    ),
    thermoPtr_
    (
        solidBlastThermo::New
        (
            fluid.mesh(),
            phaseModel::name_
        )
    ),
    rho_(thermoPtr_->rho()),
    e_(thermoPtr_->he()),
    T_(thermoPtr_->T()),
    alphaRhoPTE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTE", name_),
            fluid.mesh().time().name(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        1.5*(*this)*rho_*this->Theta_
    ),
    alphaRhoPTEPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTEPhi", name_),
            fluid.mesh().time().name(),
            fluid.mesh()
        ),
        1.5*this->alphaRhoPhi_*fvc::interpolate(Theta_)
    ),
    fluxScheme_(),
    surfTModel_(surfaceTemperatureModel::New(phaseDict_, *this))
{
    kineticTheorySystem_.addPhase(*this);
    thermoPtr_->read(phaseDict_);

    fluxScheme_ = phaseFluxScheme::NewSolid(phi_, residualAlpha().value());
    fluid.mesh().addTemporaryObject(reconstruction::ownName(alphaRho_.name()));
    fluid.mesh().addTemporaryObject(reconstruction::neiName(alphaRho_.name()));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularPhaseModel::~granularPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::ESource() const
{
    return (*this)*thermoPtr_->ESource();
}


void Foam::granularPhaseModel::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();

    //- Momentum transport
    volVectorField deltaAlphaRhoU
    (
        IOobject::groupName("deltaAlphaRhoU", name_),
        fvc::div(alphaRhoUPhi_) - alphaRho_*fluid_.g()
    );

    //- Thermal energy transport
    volScalarField deltaAlphaRhoE
    (
        IOobject::groupName("deltaAlphaRhoE", name_),
        fvc::div(alphaRhoEPhi_)
    );

    //- Pseudo thermal energy transport
    volScalarField deltaAlphaRhoPTE
    (
        IOobject::groupName("deltaAlphaRhoPTE", name_),
        fvc::div(alphaRhoPTEPhi_) + Ps_*fvc::div(phi_)
    );

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];
        if (&otherPhase != this)
        {
            if (!otherPhase.slavePressure())
            {
                deltaAlphaRhoU += (*this)*otherPhase.gradP();
            }

            if (fluid_.hasMassTransfer(*this, otherPhase))
            {
                deltaAlphaRhoU -= fluid_.mDotU(*this, otherPhase);
                deltaAlphaRhoE -= fluid_.mDotE(*this, otherPhase);
                deltaAlphaRhoPTE -= fluid_.mDotPTE(*this, otherPhase);
            }
        }
    }

    //- Solve phase mass transport
    phaseModel::solveAlphaRho();

    //- Solve thermodynamics to get energy production
    thermoPtr_->solve();
    surfTModel_->solve();

    //- Blend deltas
    deltaAlphaRhoU = cmptMultiply(deltaAlphaRhoU, solutionDs_);
    this->storeAndBlendDelta(deltaAlphaRhoU);

    //- Solve momentum transport
    this->storeAndBlendOld(alphaRhoU_);
    alphaRhoU_ -= dT*deltaAlphaRhoU;

    //- Add energy from thermodynaics
    deltaAlphaRhoE -= ESource();
    this->storeAndBlendDelta(deltaAlphaRhoE);

    // Solve thermal energy transport
    this->storeAndBlendOld(alphaRhoE_);
    alphaRhoE_ -= dT*deltaAlphaRhoE;

    //- Solve pseudo thermal energy transport
    this->storeAndBlendOld(alphaRhoPTE_);
    this->storeAndBlendDelta(deltaAlphaRhoPTE);
    alphaRhoPTE_ -= dT*(deltaAlphaRhoPTE);


    //- Update volume fraction since density is known
    alphaRho_.max(0.0);
    this->internalFieldRef() = alphaRho_()/rho_();
}


void Foam::granularPhaseModel::postExplicit()
{}


void Foam::granularPhaseModel::postImplicit()
{
    volScalarField& alpha(*this);

    if (needSolve(alpha.name()))
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

    // Solve momentum
    if (needSolve(U_.name()) || this->includeViscosity())
    {
        //- Solve momentum equation (implicit stresses)
        fvVectorMatrix UEqn
        (
            fvm::ddt(alpha, rho(), U_) - alphaRhoUAdvection_()
          + fvm::ddt(this->residualAlphaRho(), U_)
          - fvc::ddt(this->residualAlphaRho(), U_)
         ==
            models().source(alpha, rho(), U_)
        );

        if (this->includeViscosity())
        {
            // Add viscous term
            UEqn += this->divDevRhoReff(U_);
        }

        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        alphaRhoU_ = alphaRho_*U_;
    }

    // Solve thermal energy
    if (needSolve(he().name()))
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(alpha, rho(), he()) - alphaRhoEAdvection_()
          + fvm::ddt(this->residualAlphaRho(), he())
          - fvc::ddt(this->residualAlphaRho(), he())
        ==
            models().source(alpha, rho(), he())
        );
        constraints().constrain(eEqn);
        eEqn.solve();
        constraints().constrain(he());

        alphaRhoE_ = alphaRho_*he();
    }

    //- Solve granular temperature equation including solid stress and
    //  conductivity
    if (needSolve(Theta_.name()) || this->includeViscosity())
    {
        fvScalarMatrix ThetaEqn
        (
            1.5
           *(
                fvm::ddt(alpha, rho(), Theta_)
              + fvm::ddt(this->residualAlphaRho(), Theta_)
              - fvc::ddt(this->residualAlphaRho(), Theta_)
            )
          - alphaRhoPTEAdvection_()
         ==
            models().source(alpha, rho(), Theta_)
        );

        //- Solve for collisional viscosity terms
        if (this->includeViscosity())
        {
            tmp<volTensorField> tgradU
            (
                fvc::grad(fluxScheme_->interpolate(U_, U_.name()))
            );
            const volTensorField& gradU(tgradU());
            volSymmTensorField D(symm(gradU));

            volSymmTensorField tau
            (
                rho_
               *(
                    2.0*this->nut_*D
                  + (this->lambda_ - (2.0/3.0)*this->nut_)*tr(D)*I
                )
            );

            // Add solid stress and conductivity
            ThetaEqn -=
                fvm::laplacian
                (
                    this->kappa_,
                    Theta_,
                    "laplacian(kappa,Theta)"
                )
              + ((tau*alpha) && gradU);
        }
        constraints().constrain(ThetaEqn);
        ThetaEqn.solve();
        constraints().constrain(Theta_);
        Theta_.max(0);

        alphaRhoPTE_ = 1.5*alphaRho_*Theta_;
    }

    thermoPtr_->postImplicit();
    surfTModel_->postImplicit();
    dPtr_->postImplicit();
}


void Foam::granularPhaseModel::update()
{
    fluxScheme_->update
    (
        *this,
        rho_,
        U_,
        e_,
        Ptot_,
        speedOfSound()(),
        phi_,
        alphaRhoPhi_,
        alphaRhoUPhi_,
        alphaRhoEPhi_
    );

    //- Calculate PTE flux by using Riemann flux scheme to interpolate
    //  granular energy
    alphaRhoPTEPhi_ =
        1.5*fluxScheme_->flux(Theta_, alphaRho_, phi_, false);

    thermoPtr_->update();
    surfTModel_->update();
    phaseModel::update();
}


void Foam::granularPhaseModel::scaleVolumeFraction
(
    const scalar sumAlpha,
    const label celli
)
{
    (*this)[celli] /= sumAlpha;
    alphaRho_[celli] /= sumAlpha;
}


void Foam::granularPhaseModel::correctVolumeFraction
(
    const scalar alpha,
    const label celli
)
{
    (*this)[celli] = alpha;
}


void Foam::granularPhaseModel::decode()
{
    const fvConstraints& constraints = this->constraints();
    const volScalarField& alpha = *this;

    //- Correct phase mass at boundaries
    alphaRho_.correctBoundaryConditions();
    alphaRho_.boundaryFieldRef() ==
        alpha.boundaryField()*rho_.boundaryField();

    //- Store limited phase mass (only used for division)
    volScalarField alphaRhoLimited(Foam::max(alpha, residualAlpha())*rho_);

    //- Calculate velocity from momentum
    U_.internalFieldRef() = alphaRhoU_()/alphaRhoLimited();
    if (constraints.constrainsField(U_.name()))
    {
        constraints.constrain(U_);
        alphaRhoU_.internalFieldRef() = (*this)()*rho_()*U_;
    }
    U_.correctBoundaryConditions();

    //- Correct momentum at boundaries
    alphaRhoU_.correctBoundaryConditions();
    alphaRhoU_.boundaryFieldRef() ==
        alphaRho_.boundaryField()*U_.boundaryField();

    //- Limit and update thermal energy
    alphaRhoE_.max(0.0);
    e_.internalFieldRef() = alphaRhoE_()/alphaRhoLimited();
    constraints.constrain(e_);
    e_.correctBoundaryConditions();

    //- Compute granular temperature
    alphaRhoPTE_.max(0.0);
    Theta_.internalFieldRef() = alphaRhoPTE_()/(1.5*alphaRhoLimited());
    if (constraints.constrainsField(Theta_.name()))
    {
        constraints.constrain(Theta_);
        alphaRhoPTE_.internalFieldRef() = 1.5*alpha()*rho_()*Theta_();
    }
    Theta_.correctBoundaryConditions();

    alphaRhoPTE_.correctBoundaryConditions();
    alphaRhoPTE_.boundaryFieldRef() ==
        1.5*Theta_.boundaryField()*alphaRho_.boundaryField();

    thermoPtr_->correct();

    // Update total energy because e may have changed
    alphaRhoE_ == alphaRho_*e_;

    kineticTheoryModel::correct();
}


void Foam::granularPhaseModel::encode()
{
    alphaRho_ = (*this)*rho_;
    alphaRhoU_ = alphaRho_*U_;
    alphaRhoE_ = alphaRho_*e_;
    alphaRhoPTE_ = 1.5*alphaRho_*Theta_;
}


void Foam::granularPhaseModel::storeExplicit()
{
    const volScalarField& alpha = *this;
    if (needSolve(alpha.name()))
    {
        alphaAdvection_ = fvc::ddt(alpha);
    }

    if (needSolve(rho().name()))
    {
        alphaRhoAdvection_ = fvc::ddt(alphaRho_);
    }

    // Solve momentum
    if (needSolve(U_.name()) || this->includeViscosity())
    {
        alphaRhoUAdvection_ = fvc::ddt(alphaRhoU_);
    }

    // Solve thermal energy
    if (needSolve(he().name()))
    {
        alphaRhoEAdvection_ = fvc::ddt(alphaRhoE_);
    }

    //- Solve granular temperature equation including solid stress and
    //  conductivity
    if (needSolve(Theta_.name()) || this->includeViscosity())
    {
        alphaRhoPTEAdvection_ = fvc::ddt(alphaRhoPTE_);
    }

    thermoPtr_->storeExplicit();
}


void Foam::granularPhaseModel::postUpdate()
{
    thermoPtr_->postUpdate();
}


void Foam::granularPhaseModel::clear()
{
    fluxScheme_->clear();
    alphaRhoAdvection_.clear();
    alphaRhoUAdvection_.clear();
    alphaRhoEAdvection_.clear();
    alphaRhoPTEAdvection_.clear();
    thermoPtr_->clear();
}


Foam::tmp<Foam::volVectorField>
Foam::granularPhaseModel::gradP() const
{
    return fvc::grad(fluxScheme_->pf()());
}


Foam::tmp<Foam::volVectorField>
Foam::granularPhaseModel::gradAlpha() const
{
    return fvc::grad(fluxScheme_->alphaf()());
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::dissipationSource
(
    const phaseModel& phase2,
    const dimensionedScalar& deltaT
) const
{
    return
        kineticTheorySystem_.dissipationSource(*this, phase2, deltaT)
      + this->cohesion_->dissipationSource(deltaT);
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::productionSource(const phaseModel& phase) const
{
    return kineticTheorySystem_.productionSource(*this, phase);
}



Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::speedOfSound() const
{
    // Speed of sound based on particle collisions
    tmp<volScalarField> cSqr
    (
        this->pPrime()/rho_
      + 2.0/3.0
       *sqr
        (
            kineticTheorySystem_.dPsdTheta(*this)
          + cohesion_->dPsdTheta()
        )*Theta_
       /sqr(Foam::max(*this, residualAlpha())*rho_)
    );
    cSqr.ref().max(small);
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("speedOfSound", name()),
            sqrt(cSqr)
        )
    );
}

void Foam::granularPhaseModel::correctThermo()
{
    thermoPtr_->correct();
    kineticTheoryModel::correct();
}


// ************************************************************************* //
