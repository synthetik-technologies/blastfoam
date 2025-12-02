/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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

#include "coupledSinglePhaseCompressibleSystem.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledSinglePhaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        coupledSinglePhaseCompressibleSystem,
        coupled
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledSinglePhaseCompressibleSystem::coupledSinglePhaseCompressibleSystem
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coupledCompressibleSystem(dict, mesh),
    singlePhaseCompressibleSystem(dict, mesh, false)
{
    thermoPtr_->initializeModels();
    this->setModels();
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledSinglePhaseCompressibleSystem::~coupledSinglePhaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledSinglePhaseCompressibleSystem::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();

    volScalarField deltaRho
    (
        "deltaRho",
        fvc::div(rhoPhi_) // alphaRhoPhi_
    );
    this->fvTimeInt_->addDeltaSource(rho_.name(), deltaRho);

    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) // alphaRhoUPhi_
      - p_*fvc::grad(fluxScheme_->interpolate(volumeFraction_, "alpha"))
      - g_*alphaRho_
    );
    this->fvTimeInt_->addDeltaSource(rhoU_.name(), deltaRhoU);

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_) // alphaRhoEPhi
      - volumeFraction_*(rhoU_ & g_)
    );
    this->fvTimeInt_->addDeltaSource(rhoE_.name(), deltaRhoE);

    this->storeAndBlendOld(alphaRho_);
    this->storeAndBlendDelta(deltaRho);
    alphaRho_.storePrevIter();
    alphaRho_ -= dT*deltaRho;

    thermoPtr_->solve();
    deltaRhoE -= thermoPtr_->ESource();

    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendDelta(deltaRhoU);
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);

    this->storeAndBlendOld(rhoE_);
    this->storeAndBlendDelta(deltaRhoE);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::coupledSinglePhaseCompressibleSystem::postImplicit()
{
    singlePhaseCompressibleSystem::postImplicit();
}


void Foam::coupledSinglePhaseCompressibleSystem::decode()
{
    if (alphadPtr_.valid())
    {
        volumeFraction_ = min(1.0, max(0.0, 1.0 - alphadPtr_()));
        volumeFraction_.correctBoundaryConditions();
    }

    // Update density
    rho_.internalFieldRef() = alphaRho_()/max(volumeFraction_(), 1e-10);
    rho_.correctBoundaryConditions();
    alphaRho_.boundaryFieldRef() ==
        rho_.boundaryField()*volumeFraction_.boundaryField();

    // Update velocity
    volScalarField alphaRhos(alphaRho_);
    alphaRhos.max(1e-10);
    U_.internalFieldRef() = rhoU_()/alphaRhos();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() =
        alphaRho_.boundaryField()*U_.boundaryField();

    //- Update internal energy
    e_.internalFieldRef() = rhoE_()/alphaRhos() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();

    thermoPtr_->correct();

    // Update total energy since e may have changed
    rhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}


void Foam::coupledSinglePhaseCompressibleSystem::encode()
{
    if (alphadPtr_.valid())
    {
        volumeFraction_ = min(1.0, max(0.0, 1.0 - alphadPtr_()));
        volumeFraction_.correctBoundaryConditions();
    }

    alphaRho_ = rho_*volumeFraction_;
    rhoU_ = alphaRho_*U_;
    rhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}

// ************************************************************************* //
