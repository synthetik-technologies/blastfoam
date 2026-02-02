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

#include "compressibleBlastSystem.H"
#include "uniformDimensionedFields.H"
#include "fvm.H"
#include "wedgeFvPatch.H"
#include "blastRadiationModel.H"

// * * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * //

void Foam::compressibleBlastSystem::setModels()
{
    compressibleSystem::setModels();

    typeIOobject<IOdictionary> radPropertiesIO
    (
        "radiationProperties",
        rho_.time().constant(),
        rho_.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );
    if (radPropertiesIO.headerOk())
    {
        radiation_.set(blastRadiationModel::New(this->T()).ptr());
    }
}


void Foam::compressibleBlastSystem::addSources
(
    volVectorField::Internal& rhoUSource,
    volScalarField::Internal& rhoESource
) const
{

    compressibleSystem::addSources(rhoUSource, rhoESource);

    rhoESource -= thermoPtr_->ESource();
}



void Foam::compressibleBlastSystem::addESource(fvScalarMatrix& EEqn) const
{
    compressibleSystem::addESource(EEqn);

    if (radiation_.valid())
    {
        EEqn += radiation_->Sh(thermo(), e_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleBlastSystem::compressibleBlastSystem
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& thermoType
)
:
    compressibleSystem(dict, mesh),
    thermoPtr_
    (
        fluidBlastThermo::New(mesh, dict, thermoType)
    ),
    rho_(thermoPtr_().rhoRef()),
    p_(thermoPtr_->p()),
    T_(thermoPtr_->T()),
    e_(thermoPtr_->he())
{
    thermoPtr_->validate("compressibleBlastSystem", "e");
}


Foam::compressibleBlastSystem::compressibleBlastSystem
(
    const dictionary& dict,
    const fvMesh& mesh,
    autoPtr<fluidBlastThermo> thermoPtr
)
:
    compressibleSystem(dict, mesh),
    thermoPtr_(thermoPtr),
    rho_(thermoPtr_().rhoRef()),
    p_(thermoPtr_->p()),
    T_(thermoPtr_->T()),
    e_(thermoPtr_->he())
{
    thermoPtr_->validate("compressibleBlastSystem", "e");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleBlastSystem::~compressibleBlastSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::compressibleBlastSystem::update()
{
    compressibleSystem::update();
    thermoPtr_->update();
}


void Foam::compressibleBlastSystem::decode()
{
    U_.internalFieldRef() = rhoU_()/rhoEff()();
    U_.correctBoundaryConditions();

    K_ = 0.5*magSqr(U_);

    rhoU_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() =
        rhoEff().boundaryField()*U_.boundaryField();

    e_.internalFieldRef() = rhoE_()/rhoEff()() - K_();
    e_.correctBoundaryConditions();
    thermoPtr_->correct();

    //- Update total energy because the e field may have been modified
    rhoE_ = rhoEff()*(e_ + K_);
}


void Foam::compressibleBlastSystem::solve()
{
    // Update mass fields (i.e. density and volume fraction)
    solveMass();


    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU("deltaRhoU", fvc::div(rhoUPhi_));
    this->fvTimeInt_->addDeltaSource(rhoU_.name(), deltaRhoU);

    volScalarField deltaRhoE("deltaRhoE", fvc::div(rhoEPhi_));
    this->fvTimeInt_->addDeltaSource(rhoE_.name(), deltaRhoE);

    this->addSources(deltaRhoU, deltaRhoE);

    //- Store old values
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    if (this->LTS())
    {
        deltaRhoU /= corDeltaT();
        deltaRhoE /= corDeltaT();
    }

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);


    //- Solve for momentum and energy
    dimensionedScalar dT = rho_.time().deltaT();
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);
    rhoE_ -= dT*deltaRhoE;

    // Solve thermo
    thermoPtr_->solve();
}


void Foam::compressibleBlastSystem::solveImplicit()
{
    if (radiation_.valid())
    {
        radiation_->correct();
    }
    compressibleSystem::solveImplicit();

    this->thermo().solveImplicit();
    this->thermo().correct();
    constraints().constrain(p_);
    p_.correctBoundaryConditions();
}


void Foam::compressibleBlastSystem::storeExplicit()
{
    compressibleSystem::storeExplicit();
    thermoPtr_->storeExplicit();
}


void Foam::compressibleBlastSystem::clear()
{
    compressibleSystem::clear();
    thermoPtr_->clear();
}

// ************************************************************************* //
