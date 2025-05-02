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

    rhoU_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() =
        rhoEff().boundaryField()*U_.boundaryField();

    e_.internalFieldRef() = rhoE_()/rhoEff()() - 0.5*magSqr(U_());
    thermoPtr_->correct();

    //- Update total energy because the e field may have been modified
    rhoE_ = rhoEff()*(e_ + 0.5*magSqr(U_));
}


void Foam::compressibleBlastSystem::solve()
{
    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU("deltaRhoU", fvc::div(rhoUPhi_));
    this->fvTimeInt_->addDeltaSource(rhoU_.name(), deltaRhoU);

    volScalarField deltaRhoE("deltaRhoE", fvc::div(rhoEPhi_));
    this->fvTimeInt_->addDeltaSource(rhoE_.name(), deltaRhoE);

    this->addSources(deltaRhoU, deltaRhoE);

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
    if (needSolve(U_.name()) || turbulence_.valid() || dragSource_.valid())
    {
        // Solve momentum
        fvVectorMatrix UEqn
        (
            fvm::ddt(rhoEff(), U_) - fvc::ddt(rhoU_)
         ==
            models().source(rhoEff(), U_)
        );

        if (dragSource_.valid())
        {
            UEqn -= dragSource_;
        }

        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
        }
        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        rhoU_ = rhoEff()*U_;
    }

    // Solve thermal energy diffusion
    if
    (
        needSolve(e_.name())
     || turbulence_.valid()
     || radiation_.valid()
     || extESource_.valid()
    )
    {
//         if (radiation_.valid())
//         {
//             radiation_->correct();
//             rhoE_ =
//                 radiation_->calcRhoE
//                 (
//                     rho_.mesh().time().deltaT(),
//                     rhoE_,
//                     rhoEff(),
//                     e_,
//                     this->thermo().Cv()
//                 );
//         }

        if (dragSource_.valid())
        {
            rhoE_ += (dragSource_ & U_) & U_;
        }

        if (turbulence_.valid())
        {
            rhoE_ -=
                rho_.mesh().time().deltaT()
               *(U_ & fvc::div(turbulence_->devTau()));
               // *fvc::div
               //  (
               //      fvc::dotInterpolate(rho_.mesh().Sf(), turbulence_->devTau())
               //    & fluxScheme_->Uf()
               //  );
        }
        e_ = rhoE_/rhoEff() - 0.5*magSqr(U_);

        fvScalarMatrix eEqn
        (
            fvm::ddt(rhoEff(), e_) - fvc::ddt(rhoEff().prevIter(), e_)
         ==
            models().source(rhoEff(), e_)
        );
        if (extESource_.valid())
        {
            eEqn -= extESource_;
        }
        if (thermophysicalTransport_.valid())
        {
            eEqn += thermophysicalTransport_->divq(e_);
        }
        if (radiation_.valid())
        {
            radiation_->correct();
            eEqn += radiation_->Sh(thermo(), e_);
        }
        constraints().constrain(eEqn);
        eEqn.solve();
        constraints().constrain(e_);

        rhoE_ = rhoEff()*(e_ + 0.5*magSqr(U_));
    }

    this->thermo().postUpdate();
    this->thermo().correct();
    constraints().constrain(p_);
    p_.correctBoundaryConditions();

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->correct();
    }
}

// ************************************************************************* //
