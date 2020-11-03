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
#include "phaseCompressibleTurbulenceModel.H"
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
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid_.p(),
        fluid_.p().boundaryField()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid_.mesh()
    ),
    thermo_
    (
        fluidThermoModel::New
        (
            phaseName,
            p_,
            rho_,
            e_,
            T_,
            phaseDict_,
            true,
            phaseName
        )
    ),
    fluxScheme_(fluxScheme::New(fluid.mesh(), name_))
{
    thermo_->read(phaseDict_);
    thermo_->correct();

    this->turbulence_ =
        phaseCompressibleTurbulenceModel::New
        (
            *this,
            rho_,
            U_,
            alphaRhoPhi_,
            phi_,
            *this
        );
    phaseModel::initializeModels();
    thermo_->initializeModels();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::~fluidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidPhaseModel::ESource() const
{
    return (*this)*thermo_->ESource();
}


void Foam::fluidPhaseModel::solve()
{
    // Solve phase density transport to store old time values
    phaseModel::solveAlphaRho();

    // Solve thermodynamic models (activation and afterburn)
    thermo_->solve();

    // Solve momentum and energy transport
    phaseModel::solve();
}


void Foam::fluidPhaseModel::postUpdate()
{
    phaseModel::postUpdate();
    thermo_->postUpdate();
}


void Foam::fluidPhaseModel::clearODEFields()
{
    phaseModel::clearODEFields();
    thermo_->clearODEFields();
    fluxScheme_->clear();
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
            thermo_->speedOfSound()(),
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
            thermo_->speedOfSound()(),
            phi_,
            alphaRhoPhi_,
            alphaRhoUPhi_,
            alphaRhoEPhi_
        );
    }
}


void Foam::fluidPhaseModel::decode()
{
    this->correctBoundaryConditions();
    volScalarField alpha(Foam::max(*this, residualAlpha()));

    rho_.ref() = alphaRho_()/alpha();
    rho_.correctBoundaryConditions();

    alphaRho_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField();
    volScalarField alphaRho(alpha*rho_);
    alphaRho.max(1e-10);

    U_.ref() = alphaRhoU_()/(alphaRho());
    U_.correctBoundaryConditions();

    alphaRhoU_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField()*U_.boundaryField();

    volScalarField E(alphaRhoE_/alphaRho);
    e_.ref() = E() - 0.5*magSqr(U_());

    //--- Hard limit, e
    if(Foam::min(e_).value() < 0)
    {
        e_.max(small);
        alphaRhoE_.ref() = (*this)()*rho_()*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    alphaRhoE_.boundaryFieldRef() =
        (*this).boundaryField()
       *rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_->correct();
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
    thermo_->correct();
}

// ************************************************************************* //
