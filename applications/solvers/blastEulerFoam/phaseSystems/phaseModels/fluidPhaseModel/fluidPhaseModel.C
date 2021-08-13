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
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
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
    const label index
)
:
    phaseModel(fluid, phaseName, index),
    thermoPtr_
    (
        fluidBlastThermo::New
        (
            1,
            fluid.mesh(),
            phaseDict_,
            phaseName
        )
    ),
    rho_(thermoPtr_->rho()),
    e_(thermoPtr_->he()),
    T_(thermoPtr_->T()),
    p_(thermoPtr_->p()),
    fluxScheme_(phaseFluxScheme::New(fluid.mesh(), name_))
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

    phaseModel::initializeModels();
    thermoPtr_->initializeModels();
    thermoPtr_->correct();
}


Foam::fluidPhaseModel::fluidPhaseModel
(
    const label nPhases,
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    phaseModel(fluid, phaseName, index),
    thermoPtr_
    (
        fluidBlastThermo::New
        (
            nPhases,
            fluid.mesh(),
            phaseDict_,
            phaseName
        )
    ),
    rho_(thermoPtr_->rho()),
    e_(thermoPtr_->he()),
    T_(thermoPtr_->T()),
    p_(thermoPtr_->p()),
    fluxScheme_(phaseFluxScheme::New(fluid.mesh(), name_))
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
    // Solve phase density transport to store old time values
    phaseModel::solveAlphaRho();

    // Solve thermodynamic models (activation and afterburn)
    thermoPtr_->solve();

    // Solve momentum and energy transport
    phaseModel::solve();
}


void Foam::fluidPhaseModel::postUpdate()
{
    phaseModel::postUpdate();
    thermoPtr_->postUpdate();
}


void Foam::fluidPhaseModel::update()
{
    const volScalarField& alpha = *this;
    if (this->solveAlpha_)
    {
        fluxScheme_->update
        (
            alpha,
            rho_,
            U_,
            e_,
            p_,
            thermoPtr_->speedOfSound()(),
            phi_,
            this->alphaPhiPtr_(),
            alphaRhoPhi_,
            alphaRhoUPhi_,
            alphaRhoEPhi_
        );
    }
    else
    {
        fluxScheme_->update
        (
            alpha,
            rho_,
            U_,
            e_,
            p_,
            thermoPtr_->speedOfSound()(),
            phi_,
            alphaRhoPhi_,
            alphaRhoUPhi_,
            alphaRhoEPhi_
        );
    }

    thermoPtr_->update();
    phaseModel::update();
}


void Foam::fluidPhaseModel::decode()
{
    const fvConstraints& constraints(this->fluid_.constraints());

    this->correctBoundaryConditions();
    if (constraints.constrainsField(this->name()))
    {
        constraints.constrain(*this);
    }
    volScalarField alpha(Foam::max(*this, residualAlpha()));


    rho_.ref() = alphaRho_()/alpha();
    if (constraints.constrainsField(rho_.name()))
    {
        constraints.constrain(rho_);
        alphaRho_.ref() = (*this)()*rho_;
    }
    rho_.correctBoundaryConditions();
    alphaRho_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField();
    volScalarField alphaRho(alpha*rho_);
    alphaRho.max(1e-10);

    U_.ref() = alphaRhoU_()/(alphaRho());
    if (constraints.constrainsField(U_.name()))
    {
        constraints.constrain(U_);
        alphaRhoU_.ref() = alphaRho_()*U_();
    }
    U_.correctBoundaryConditions();

    alphaRhoU_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField()*U_.boundaryField();

    volScalarField E(alphaRhoE_/alphaRho);
    e_.ref() = E() - 0.5*magSqr(U_());
    if (constraints.constrainsField(e_.name()))
    {
        constraints.constrain(e_);
        alphaRhoE_.ref() = alphaRho_()*(e_() + 0.5*magSqr(U_()));
    }

    //- Limit internal energy it there is a negative temperature
    if (Foam::min(this->T_).value() < 0.0 && thermoPtr_->limit())
    {
        if (debug)
        {
            WarningInFunction
                << "Lower limit of temperature reached, min(T) = "
                << Foam::min(T_).value()
                << ", limiting internal energy." << endl;
        }
        volScalarField limit(pos(T_));
        T_.max(small);
        e_ = e_*limit + thermoPtr_->he(p_, T_)*(1.0 - limit);
        alphaRhoE_.ref() = alphaRho_()*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    alphaRhoE_.boundaryFieldRef() =
        (*this).boundaryField()
       *rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );
    thermoPtr_->correct();

    if (constraints.constrainsField(p_.name()))
    {
        constraints.constrain(p_);
        p_.correctBoundaryConditions();
    }
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
