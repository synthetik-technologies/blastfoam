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
    thermo_
    (
        solidThermoModel::New
        (
            phaseName,
            fluid.mesh(),
            phaseDict_,
            true,
            phaseName
        )
    ),
    rho_(thermo_->rho()),
    p_(fluid.p()),
    alphaRhoPTE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTE", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.5*(*this)*rho_*this->Theta_
    ),
    alphaRhoPTEPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTEPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.5*this->alphaRhoPhi_*fvc::interpolate(Theta_)
    ),
    residualAlpha_
    (
        thermo_->residualAlpha()
    ),
    fluxScheme_(new fluxSchemes::AUSMPlusUp(fluid.mesh(), phaseName))
{
    thermo_->read(phaseDict_);
    phaseModel::initializeModels();
    thermo_->initializeModels();
    encode();
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
        fvc::div(alphaRhoUPhi_) - alphaRho_*fluid_.g()
    );
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];
        if (!phase.granular())
        {
            deltaAlphaRhoU += (*this)*phase.gradP();
        }
    }

    //- Thermal energy transport
    volScalarField deltaAlphaRhoE(fvc::div(alphaRhoEPhi_));

    //- Pseudo thermal energy transport
    volScalarField deltaAlphaRhoPTE(fvc::div(alphaRhoPTEPhi_) + Ps_*fvc::div(phi_));

    //- Solve phase mass transport
    phaseModel::solveAlphaRho();

    //- Solve thermodynamics to get energy production
    thermo_->solve();

    //- Blend deltas
    this->storeAndBlendDelta(deltaAlphaRhoU, deltaAlphaRhoU_);

    //- Add energy from thermodynaics
    deltaAlphaRhoE -= ESource();
    this->storeAndBlendDelta(deltaAlphaRhoE, deltaAlphaRhoE_);

    this->storeAndBlendDelta(deltaAlphaRhoPTE, deltaAlphaRhoPTE_);

    //- Solve momentum transport
    this->storeAndBlendOld(alphaRhoU_, alphaRhoUOld_);
    alphaRhoU_ -= cmptMultiply(dT*deltaAlphaRhoU, solutionDs_);
    alphaRhoU_.correctBoundaryConditions();

    // Solve thermal energy transport
    this->storeAndBlendOld(alphaRhoE_, alphaRhoEOld_);
    alphaRhoE_ -= dT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();

    //- Solve pseudo thermal energy transport
    this->storeAndBlendOld(alphaRhoPTE_, alphaRhoPTEOld_);
    alphaRhoPTE_ -= dT*(deltaAlphaRhoPTE);
    alphaRhoPTE_.correctBoundaryConditions();

    //- Update diameterModel
    this->dPtr_->solve();
}


void Foam::granularPhaseModel::postUpdate()
{
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
        dimensionedScalar smallAlphaRho(dimDensity, 1e-10);

        //- Solve momentum equation (implicit stresses)
        fvVectorMatrix UEqn
        (
            fvm::ddt(alphaRho_, U_) - fvc::ddt(alphaRho_, U_)
          + fvm::ddt(smallAlphaRho, U_) - fvc::ddt(smallAlphaRho, U_)
          + this->divDevRhoReff(U_)
        );
        UEqn.solve();

        //- Solve granular temperature equation including solid stress and
        //  conductivity
        fvScalarMatrix ThetaEqn
        (
            1.5
           *(
                fvm::ddt(alphaRho_, Theta_) - fvc::ddt(alphaRho_, Theta_)
              + fvm::ddt(smallAlphaRho, Theta_) - fvc::ddt(smallAlphaRho, Theta_)
            )
          - fvm::laplacian(this->kappa_, Theta_)
         ==
            ((tau*(*this)) && gradU)
        );
        ThetaEqn.solve();

        //- Update conservative variables
        alphaRhoU_ =  cmptMultiply(alphaRho_*U_, solutionDs_);
        alphaRhoPTE_ =  alphaRho_*1.5*Theta_;

    }

    dPtr_->postUpdate();
    thermo_->postUpdate();
}


void Foam::granularPhaseModel::clearODEFields()
{
    phaseModel::clearODEFields();
    fluxScheme_->clear();

    this->clearOld(alphaRhoPTEOld_);
    this->clearDelta(deltaAlphaRhoPTE_);

    thermo_->clearODEFields();
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

    thermo_->update();
    phaseModel::update();
}


void Foam::granularPhaseModel::decode()
{
    volScalarField& alpha(*this);

    //- Update volume fraction since density is known
    alpha.ref() = alphaRho_()/rho();
    alpha.correctBoundaryConditions();

    //- Correct phase mass at boundaries
    alphaRho_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField();

    //- Store limited phase mass (only used for division)
    volScalarField alphaRho(Foam::max(alpha, residualAlpha_)*rho_);

    //- Calculate velocity from momentum
    U_.ref() = alphaRhoU_()/alphaRho();
    U_.correctBoundaryConditions();

    //- Correct momentum at boundaries
    alphaRhoU_.boundaryFieldRef() = alphaRho_.boundaryField()*U_.boundaryField();

    //- Limit and update thermal energy
    alphaRhoE_.max(0.0);
    e_.ref() = alphaRhoE_()/alphaRho();
    e_.correctBoundaryConditions();
    alphaRhoE_.boundaryFieldRef() = alphaRho_.boundaryField()*e_.boundaryField();

    //- Compute granular temperature
    alphaRhoPTE_.max(0.0);
    Theta_.ref() = alphaRhoPTE_()/(1.5*alphaRho());
    Theta_.correctBoundaryConditions();
    alphaRhoPTE_.boundaryFieldRef() =
        1.5*Theta_.boundaryField()*alphaRho_.boundaryField();

    thermo_->correct();
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
    T_ = thermo_->calcT();
    kineticTheoryModel::correct();
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::ESource() const
{
    return (*this)*thermo_->ESource();
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::alpha() const
{
    return thermo_->alpha();
}


Foam::tmp<Foam::scalarField>
Foam::granularPhaseModel::alpha(const label patchi) const
{
    return thermo_->alpha(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::alphahe() const
{
    return thermo_->alphahe();
}


Foam::tmp<Foam::scalarField>
Foam::granularPhaseModel::alphahe(const label patchi) const
{
    return thermo_->alphahe(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::Cp() const
{
    return thermo_->Cp();
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::Cv() const
{
    return thermo_->Cv();
}

// ************************************************************************* //
