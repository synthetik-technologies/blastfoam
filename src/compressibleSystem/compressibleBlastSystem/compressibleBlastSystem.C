/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
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
    U_.ref() = rhoU_()/rho()();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() =
        rho().boundaryField()*U_.boundaryField();

    e_.ref() = rhoE_()/rho()() - 0.5*magSqr(U_());

    thermoPtr_->correct();

    //- Update total energy because the e field may have been modified
    rhoE_ = rho()*(e_ + 0.5*magSqr(U_));
}


void Foam::compressibleBlastSystem::solve()
{
    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) - g_*rho()
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
    bool updateE = false;
    if (radiation_.valid())
    {
        updateE = true;
        radiation_->correct();
        rhoE_ =
            radiation_->calcRhoE
            (
                rho_.mesh().time().deltaT(),
                rhoE_,
                rho(),
                e_,
                this->thermo().Cv()
            );
    }

    if (needSolve(U_.name()) || turbulence_.valid())
    {
        updateE = true;

        // Solve momentum
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho(), U_) - fvc::ddt(rhoU_)
         ==
            models().source(rho(), U_)
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
        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        rhoU_ = rho()*U_;
    }

    //- Update internal energy
    if (updateE)
    {
        e_ = rhoE_/rho() - 0.5*magSqr(U_);
        e_.correctBoundaryConditions();
    }

    // Solve thermal energy diffusion
    if (needSolve(e_.name()) || turbulence_.valid())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho(), e_) - fvc::ddt(rho().prevIter(), e_)
         ==
            models().source(rho(), e_)
        );
        if (extESource_.valid())
        {
            eEqn -= extESource_;
        }
        if (turbulence_.valid())
        {
            eEqn += thermophysicalTransport_->divq(e_);
        }
        constraints().constrain(eEqn);
        eEqn.solve();
        constraints().constrain(e_);

        rhoE_ = rho()*(e_ + 0.5*magSqr(U_));
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    this->thermo().postUpdate();
    this->thermo().correct();
    constraints().constrain(p_);
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
