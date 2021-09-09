/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

void Foam::compressibleBlastSystem::setModels()
{
    compressibleSystem::setModels();

    IOobject radPropertiesIO
    (
        "radiationProperties",
        rho_.time().constant(),
        rho_.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );
    if (radPropertiesIO.typeHeaderOk<IOdictionary>(true))
    {
        needPostUpdate_ = true;
        radiation_.set(blastRadiationModel::New(this->T()).ptr());
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleBlastSystem::compressibleBlastSystem
(
    const label nPhases,
    const fvMesh& mesh
)
:
    compressibleSystem(mesh),
    thermoPtr_
    (
        fluidBlastThermo::New(nPhases, mesh, *this)
    ),
    rho_(thermoPtr_->rho()),
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
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    e_.ref() = rhoE_()/rho_() - 0.5*magSqr(U_());

    thermoPtr_->correct();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );
}


void Foam::compressibleBlastSystem::solve()
{
    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) - g_*rho_
    );

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_)
      - ESource()
      - (rhoU_ & g_)
    );

    //- Store old values
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    //- Solve for momentum and energy
    dimensionedScalar dT = rho_.time().deltaT();
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::compressibleBlastSystem::postUpdate()
{
    this->thermo().postUpdate();

    if (!needPostUpdate_)
    {
        return;
    }

    if (radiation_.valid())
    {
        radiation_->correct();
        rhoE_ =
            radiation_->calcRhoE
            (
                rho_.mesh().time().deltaT(),
                rhoE_,
                rho_,
                e_,
                this->thermo().Cv()
            );
    }

    if (solveFields_.found(U_.name()) || turbulence_.valid())
    {
        // Solve momentum
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
         ==
            modelsPtr_->source(rho_, U_)
        );

        if (dragSource_.valid())
        {
            UEqn -= dragSource_;
        }

        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
            rhoE_ +=
                rho_.mesh().time().deltaT()
               *fvc::div
                (
                    fvc::dotInterpolate(rho_.mesh().Sf(), turbulence_->devTau())
                  & fluxScheme_->Uf()
                );
        }
        constraintsPtr_->constrain(UEqn);
        UEqn.solve();
        constraintsPtr_->constrain(U_);

        // Update internal energy
        e_ = rhoE_/rho_ - 0.5*magSqr(U_);
    }

    // Solve thermal energy diffusion
    if (solveFields_.found(e_.name()) || turbulence_.valid())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho_, e_) - fvc::ddt(rho_, e_)
         ==
            modelsPtr_->source(rho_, e_)
        );
        if (extESource_.valid())
        {
            eEqn -= extESource_;
        }
        if (turbulence_.valid())
        {
            eEqn += thermophysicalTransport_->divq(e_);
        }
        constraintsPtr_->constrain(eEqn);
        eEqn.solve();
        constraintsPtr_->constrain(e_);
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    encode();
    this->thermo().correct();
    constraintsPtr_->constrain(p_);
}


void Foam::compressibleBlastSystem::addECoeff
(
    const volScalarField::Internal& ECoeff
)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(e_, dimEnergy/dimTime)
            );
    }
    extESource_.ref() -= fvm::Sp(ECoeff, e_);
}


void Foam::compressibleBlastSystem::addESource
(
    const volScalarField::Internal& ESource
)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(e_, dimEnergy/dimTime)
            );
    }
    extESource_.ref() += ESource;
}


void Foam::compressibleBlastSystem::addUCoeff
(
    const volScalarField::Internal& UCoeff
)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() -= fvm::Sp(UCoeff, U_);
}


void Foam::compressibleBlastSystem::addUSource
(
    const volVectorField::Internal& USource
)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() += USource;
}

// ************************************************************************* //
