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

    K_ = 0.5*magSqr(U_);

    rhoU_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() =
        rhoEff().boundaryField()*U_.boundaryField();

    e_.internalFieldRef() = rhoE_()/rhoEff()() - K_();
    e_.correctBoundaryConditions();
    thermoPtr_->correct();

    //- Update total energy because the e field may have been modified
    rhoE_ = rhoEff()*(e_ + K_);

    if (explicitViscosity_)
    {
        if (turbulence_.valid())
        {
            turbulence_->predict();
        }
        if (thermophysicalTransport_.valid())
        {
            thermophysicalTransport_->predict();
        }
    }
}


void Foam::compressibleBlastSystem::solve()
{
    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU("deltaRhoU", fvc::div(rhoUPhi_));
    this->fvTimeInt_->addDeltaSource(rhoU_.name(), deltaRhoU);

    volScalarField deltaRhoE("deltaRhoE", fvc::div(rhoEPhi_));
    this->fvTimeInt_->addDeltaSource(rhoE_.name(), deltaRhoE);

    this->addSources(deltaRhoU, deltaRhoE);

    if (explicitViscosity_ && turbulence_.valid())
    {
        tmp<volSymmTensorField> tdevTau(turbulence_->devTau());
        tdevTau.ref() += (2.0/3.0)*rhoEff()*turbulence_->k()*symmTensor::I;

        deltaRhoU += fvc::div(tdevTau());
        deltaRhoE += fvc::div
            (
                fvc::dotInterpolate(mesh().Sf(), tdevTau)
              & flux().Uf()
            )
          + fvc::div(thermophysicalTransport_->q()*mesh().magSf());
    }

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
}


void Foam::compressibleBlastSystem::postUpdate()
{
    if (turbulence_.valid())
    {
        turbulence_->predict();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->predict();
    }

    tmp<surfaceVectorField> devTau;
    if
    (
        needSolve(U_.name())
     || (!explicitViscosity_ && turbulence_.valid())
     || dragSource_.valid()
    )
    {
        tmp<fvVectorMatrix> divDevTau;
        for (label iter = 0; iter < 2; iter++)
        {
            if (!explicitViscosity_ && turbulence_.valid())
            {
                divDevTau =
                    turbulence_->divDevTau(U_)
                  + fvc::grad((2.0/3.0)*rhoEff()*turbulence_->k());
            }

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
            if (divDevTau.valid())
            {
                UEqn += divDevTau();
            }

            UEqn.relax();

            constraints().constrain(UEqn);
            UEqn.solve();
            constraints().constrain(U_);
        }

        if (divDevTau.valid())
        {
            devTau = divDevTau().flux();
            // devTau =
            //     fvc::dotInterpolate
            //     (
            //         mesh().Sf(),
            //         turbulence_->devTau()
            //     )
            //   + fvc::interpolate
            //     (
            //         (2.0/3.0)*rhoEff()*turbulence_->k()
            //     )*mesh().Sf();
        }

        K_ = 0.5*magSqr(U_);
        rhoU_ = rhoEff()*U_;
    }

    // Solve thermal energy diffusion
    if
    (
        needSolve(e_.name())
     || (!explicitViscosity_ && turbulence_.valid())
     || radiation_.valid()
     || extESource_.valid()
    )
    {
        // e_ = rhoE_/rhoEff() - 0.5*magSqr(U_);
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

        fvScalarMatrix EEqn
        (
            fvm::ddt(rhoEff(), e_)
          - fvc::ddt(rhoE_) // Advection only
          + fvc::ddt(rhoEff(), K_)
         ==
            models().source(rhoEff(), e_)
        );

        if (devTau.valid())
        {
            // tmp<volScalarField> tk(turbulence_->k());
            // const volScalarField& k = tk();
            EEqn +=
                fvc::div(devTau & flux().Uf());
              // + fvc::ddt(rhoEff(), k);
        }
        if (!explicitViscosity_ && thermophysicalTransport_.valid())
        {
            EEqn += thermophysicalTransport_->divq(e_);
        }
        if (dragSource_.valid())
        {
            EEqn -= (dragSource_ & U_) & U_;
        }
        if (extESource_.valid())
        {
            EEqn -= extESource_;
        }
        if (radiation_.valid())
        {
            radiation_->correct();
            EEqn += radiation_->Sh(thermo(), e_);
        }

        EEqn.relax();

        constraints().constrain(EEqn);
        EEqn.solve();
        constraints().constrain(e_);

        rhoE_ = rhoEff()*(e_ + K_);
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
