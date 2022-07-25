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
            phaseDict_,
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
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        1.5*(*this)*rho_*this->Theta_
    ),
    alphaRhoPTEPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTEPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        1.5*this->alphaRhoPhi_*fvc::interpolate(Theta_)
    ),
    fluxScheme_(phaseFluxScheme::NewSolid(fluid.mesh(), phaseName))
{
    thermoPtr_->read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularPhaseModel::~granularPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    if (fluid_.hasMassTransfer(*this))
    {
	forAll(fluid_.phases(), phasei)
	{
	    const phaseModel& otherPhase = fluid_.phases()[phasei];
	    if (&otherPhase != this && fluid_.hasMassTransfer(*this, otherPhase))
	    {
		if (!otherPhase.slavePressure())
		{
		    deltaAlphaRhoU += (*this)*otherPhase.gradP();
		}
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

    //- Blend deltas
    deltaAlphaRhoU = cmptMultiply(deltaAlphaRhoU, solutionDs_);
    this->storeAndBlendDelta(deltaAlphaRhoU);

    //- Solve momentum transport
    this->storeAndBlendOld(alphaRhoU_);
    alphaRhoU_ -= dT*deltaAlphaRhoU;
    alphaRhoU_.correctBoundaryConditions();

    //- Add energy from thermodynaics
    deltaAlphaRhoE -= ESource();
    this->storeAndBlendDelta(deltaAlphaRhoE);

    // Solve thermal energy transport
    this->storeAndBlendOld(alphaRhoE_);
    alphaRhoE_ -= dT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();

    //- Solve pseudo thermal energy transport
    this->storeAndBlendOld(alphaRhoPTE_);
    this->storeAndBlendDelta(deltaAlphaRhoPTE);
    alphaRhoPTE_ -= dT*(deltaAlphaRhoPTE);
    alphaRhoPTE_.correctBoundaryConditions();
}


void Foam::granularPhaseModel::postUpdate()
{
    dimensionedScalar smallAlphaRho(dimDensity, 1e-6);
    volScalarField& alpha(*this);

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

    // Solve momentum
    if (needSolve(U_.name()) || this->includeViscosity())
    {
        //- Solve momentum equation (implicit stresses)
        fvVectorMatrix UEqn
        (
            fvm::ddt(alphaRho_, U_) - fvc::ddt(alphaRhoU_)
          + fvm::ddt(smallAlphaRho, U_) - fvc::ddt(smallAlphaRho, U_)
         ==
            models().source(alphaRho_, U_)
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
            fvm::ddt(alphaRho_, he()) - fvc::ddt(alphaRhoE_)
          + fvm::ddt(smallAlphaRho, he()) - fvc::ddt(smallAlphaRho, he())
        ==
            models().source(alphaRho_, he())
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
              - fvc::ddt(alphaRho_.prevIter(), Theta_)
              + fvm::ddt(smallAlphaRho, Theta_)
              - fvc::ddt(smallAlphaRho, Theta_)
            )
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

    thermoPtr_->postUpdate();
    dPtr_->postUpdate();
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
        1.5*alphaRhoPhi_*fluxScheme_->interpolate(Theta_, "Theta");

    thermoPtr_->update();
    phaseModel::update();
}


void Foam::granularPhaseModel::correctVolumeFraction()
{
    volScalarField& alpha(*this);

    //- Update volume fraction since density is known
    alphaRho_.max(0.0);
    alpha.ref() = alphaRho_()/rho_();
    alpha.correctBoundaryConditions();
}


void Foam::granularPhaseModel::decode()
{
    const volScalarField& alpha(*this);

    //- Correct phase mass at boundaries
    alphaRho_.boundaryFieldRef() ==
        alpha.boundaryField()*rho_.boundaryField();

    //- Store limited phase mass (only used for division)
    volScalarField alphaRhoLimited(Foam::max(alpha, residualAlpha())*rho_);

    //- Calculate velocity from momentum
    U_.ref() = alphaRhoU_()/alphaRhoLimited();
    U_.correctBoundaryConditions();

    //- Correct momentum at boundaries
    alphaRhoU_.boundaryFieldRef() == alphaRho_.boundaryField()*U_.boundaryField();

    //- Limit and update thermal energy
    alphaRhoE_.max(0.0);
    e_.ref() = alphaRhoE_()/alphaRhoLimited();

    //- Compute granular temperature
    alphaRhoPTE_.max(0.0);
    Theta_.ref() = alphaRhoPTE_()/(1.5*alphaRhoLimited());
    Theta_.correctBoundaryConditions();
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
Foam::granularPhaseModel::dissipationSource(const phaseModel& phase2) const
{
    // Dissipation of granular energy (Huilin and Gidaspow 2003, Eq. 25)
    scalar pi(Foam::constant::mathematical::pi);
    volScalarField Theta1(Theta_);
    Theta1.max(1e-10);
    volScalarField Theta2(phase2.Theta());
    Theta2.max(1e-10);
    phasePairKey key(name(), phase2.name(), false);

    volScalarField m1(pi/6.0*pow3(this->d())*rho_);
    volScalarField m2(pi/6.0*pow3(phase2.d())*phase2.rho());
    volScalarField m0(m1 + m2);
    volScalarField m1Thetam2Theta(sqr(m1)*Theta1 + sqr(m2)*Theta2);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            (
                3.0/this->d()
               *sqrt(2.0*sqr(m0)*Theta_*phase2.Theta()/(pi*m1Thetam2Theta))
              - (3.0*m0*(m1*Theta_ + m2*phase2.Theta()))/(4.0*m1Thetam2Theta)
               *fvc::div(phi_)
            )
           *(1.0 - kineticTheorySystem_.es(key))
           *kineticTheorySystem_.Ps(*this, phase2)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::productionSource(const phaseModel& phase) const
{
    // Production of granular energy (Houim and Oran 2016, Eq. 3.49, Eq. B 66)
    return tmp<volScalarField>
    (
        new volScalarField
        (
            81.0*(*this)*sqr(phase.mu())*magSqr(U_ - phase.U())
           /(gs0_*pow3(d())*rho_*sqrt(Foam::constant::mathematical::pi))
        )
    );
}



Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::speedOfSound() const
{
    // Speed of sound based on particle collisions
    tmp<volScalarField> cSqr
    (
        this->pPrime()/rho_
      + 2.0/3.0*sqr(kineticTheorySystem_.dPsdTheta(*this))*Theta_
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


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::ESource() const
{
    return (*this)*thermoPtr_->ESource();
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::alpha() const
{
    return thermoPtr_->alpha();
}


Foam::tmp<Foam::scalarField>
Foam::granularPhaseModel::alpha(const label patchi) const
{
    return thermoPtr_->alpha(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::alphahe() const
{
    return thermoPtr_->alphahe();
}


Foam::tmp<Foam::scalarField>
Foam::granularPhaseModel::alphahe(const label patchi) const
{
    return thermoPtr_->alphahe(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::Cp() const
{
    return thermoPtr_->Cp();
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::Cv() const
{
    return thermoPtr_->Cv();
}

// ************************************************************************* //
