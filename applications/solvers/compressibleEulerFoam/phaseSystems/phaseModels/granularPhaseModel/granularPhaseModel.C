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
            phaseDict_
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
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularPhaseModel::~granularPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::granularPhaseModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (oldIs_[stepi - 1] != -1)
    {
        alphaRhoOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(alphaRho_)
        );
        alphaRhoUOld_.set
        (
            oldIs_[stepi - 1],
            new volVectorField(alphaRhoU_)
        );
        alphaRhoEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(alphaRhoE_)
        );
        alphaRhoPTEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(alphaRhoPTE_)
        );
    }
    volScalarField alphaRhoOld(ai[stepi - 1]*(alphaRho_));
    volVectorField alphaRhoUOld(ai[stepi - 1]*(alphaRhoU_));
    volScalarField alphaRhoEOld(ai[stepi - 1]*(alphaRhoE_));
    volScalarField alphaRhoPTEOld(ai[stepi - 1]*(alphaRhoPTE_));

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            alphaRhoOld += ai[fi]*alphaRhoOld_[fi];
            alphaRhoUOld += ai[fi]*alphaRhoUOld_[fi];
            alphaRhoEOld += ai[fi]*alphaRhoEOld_[fi];
            alphaRhoPTEOld += ai[fi]*alphaRhoPTEOld_[fi];
        }
    }

    volScalarField deltaAlphaRho(fvc::div(alphaRhoPhi_));
    volVectorField deltaAlphaRhoU(fvc::div(alphaRhoUPhi_));
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];
        if (!phase.granular())
        {
            deltaAlphaRhoU += (*this)*phase.gradP();
        }
    }
    volScalarField deltaAlphaRhoE(fvc::div(alphaRhoEPhi_));
    volScalarField deltaAlphaRhoPTE
    (
        fvc::div(alphaRhoPTEPhi_)
      + Ps_*fvc::div(phi_)
    );

    if (deltaIs_[stepi - 1] != -1)
    {
        deltaAlphaRho_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRho)
        );
        deltaAlphaRhoU_.set
        (
            deltaIs_[stepi - 1],
            new volVectorField(deltaAlphaRhoU)
        );
        deltaAlphaRhoE_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRhoE)
        );
        deltaAlphaRhoPTE_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRhoPTE)
        );
    }
    deltaAlphaRho *= bi[stepi - 1];
    deltaAlphaRhoU *= bi[stepi - 1];
    deltaAlphaRhoE *= bi[stepi - 1];
    deltaAlphaRhoPTE *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaAlphaRho += bi[fi]*deltaAlphaRho_[fi];
            deltaAlphaRhoU += bi[fi]*deltaAlphaRhoU_[fi];
            deltaAlphaRhoE += bi[fi]*deltaAlphaRhoE_[fi];
            deltaAlphaRhoPTE += bi[fi]*deltaAlphaRhoPTE_[fi];
        }
    }


    dimensionedScalar dT = rho_.time().deltaT();
    vector solutionDs((vector(this->mesh().solutionD()) + vector::one)/2.0);

    alphaRho_ = alphaRhoOld - dT*(deltaAlphaRho);
    alphaRho_.max(0);
    alphaRhoU_ = cmptMultiply(alphaRhoUOld - dT*deltaAlphaRhoU, solutionDs);
    alphaRhoE_ = alphaRhoEOld - dT*(deltaAlphaRhoE);
    alphaRhoPTE_ = alphaRhoPTEOld - dT*(deltaAlphaRhoPTE);
}


void Foam::granularPhaseModel::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    phaseModel::setODEFields(nSteps, storeFields, storeDeltas);
    alphaRhoPTEOld_.resize(nOld_);
    deltaAlphaRhoPTE_.resize(nDelta_);
}

void Foam::granularPhaseModel::clearODEFields()
{
    phaseModel::clearODEFields();
    fluxScheme_->clear();
    alphaRhoPTEOld_.clear();
    deltaAlphaRhoPTE_.clear();

    alphaRhoPTEOld_.resize(nOld_);
    deltaAlphaRhoPTE_.resize(nDelta_);
}


void Foam::granularPhaseModel::update()
{
    volScalarField c(speedOfSound());
    fluxScheme_->update
    (
        *this,
        rho_,
        U_,
        e_,
        Ptot_,
        c,
        phi_,
        alphaRhoPhi_,
        alphaRhoUPhi_,
        alphaRhoEPhi_
    );
    alphaRhoPTEPhi_ =
        1.5*alphaRhoPhi_*fluxScheme_->interpolate(Theta_, "Theta");
}


void Foam::granularPhaseModel::decode()
{
    volScalarField& alpha(*this);

    alpha.ref() = alphaRho_()/rho();
    alpha.correctBoundaryConditions();
    alphaRho_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField();
    volScalarField alphaRho(Foam::max(alpha, residualAlpha_)*rho_);

    U_.ref() = alphaRhoU_()/(alphaRho());
    U_.correctBoundaryConditions();

    alphaRhoU_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField()*U_.boundaryField();

    e_.ref() = alphaRhoE_/alphaRho;

    //- Limit internal energy it there is a negative temperature
    if(Foam::min(T_).value() < small && thermo_->limit())
    {
        if (debug)
        {
            WarningInFunction
                << "Lower limit of temperature reached, min(T) = "
                << Foam::min(T_).value()
                << ", limiting internal energy." << endl;
        }
        volScalarField limit(pos(T_ - small));
        T_.max(small);
        e_ = e_*limit + thermo_->E()*(1.0 - limit);
        alphaRhoE_.ref() = alphaRho_()*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();
    alphaRhoE_.boundaryFieldRef() =
        alphaRho.boundaryField()*e_.boundaryField();

    Theta_.ref() = alphaRhoPTE_()/(1.5*alphaRho());
    Theta_.max(0.0);
    Theta_.correctBoundaryConditions();
    alphaRhoPTE_.boundaryFieldRef() =
        Theta_.boundaryField()*alphaRho.boundaryField();

    T_ = thermo_->calcT();
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
    scalar pi(Foam::constant::mathematical::pi);
    volScalarField Theta1(Theta_);
    Theta1.max(1e-10);
    volScalarField Theta2(phase2.Theta());
    Theta2.max(1e-10);
    phasePairKey key(name(), phase2.name(), false);

    volScalarField m1(pi/6.0*pow3(d_)*rho_);
    volScalarField m2(pi/6.0*pow3(phase2.d())*phase2.rho());
    volScalarField m0(m1 + m2);
    volScalarField m1Thetam2Theta(sqr(m1)*Theta1 + sqr(m2)*Theta2);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            (
                3.0/d_
               *sqrt
                (
                    2.0*sqr(m0)*Theta_*phase2.Theta()
                    /(pi*m1Thetam2Theta)
                )
              - (3.0*m0*(m1*Theta_ + m2*phase2.Theta()))
               /(4.0*m1Thetam2Theta)
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
    return tmp<volScalarField>
    (
        new volScalarField
        (
            81.0*(*this)*sqr(phase.mu())*magSqr(U_ - phase.U())
           /(
                gs0_
               *pow3(d())
               *rho_
               *sqrt(Foam::constant::mathematical::pi)
            )
        )
    );
}



Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::speedOfSound() const
{
    tmp<volScalarField> cSqr
    (
        pPrime()/rho_
      + 2.0/3.0*sqr(kineticTheorySystem_.dPsdTheta(*this))*Theta_
       /sqr(Foam::max(*this, residualAlpha())*rho_)
    );
    cSqr.ref().max(0.0);
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
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("0", alphaRhoE_.dimensions(), 0.0)
        )
    );
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
