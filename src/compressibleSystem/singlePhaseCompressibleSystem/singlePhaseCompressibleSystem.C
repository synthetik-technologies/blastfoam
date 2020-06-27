/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "singlePhaseCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        phaseCompressibleSystem,
        singlePhaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::singlePhaseCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseCompressibleSystem(mesh, dict),
    thermo_
    (
        fluidThermoModel::New
        (
            word::null,
            p_,
            rho_,
            e_,
            T_,
            dict.subDict("mixture"),
            true
        )
    )
{
    setModels(dict);
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::~singlePhaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singlePhaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    volScalarField rhoOld(rho_);
    if (rho_.mesh().moving() && stepi == 1)
    {
        rhoOld.ref() *= rho_.mesh().Vsc0()/rho_.mesh().Vsc();
    }

    if (oldIs_[stepi - 1] != -1)
    {
        rhoOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoOld)
        );
    }
    rhoOld *= ai[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoOld += ai[fi]*rhoOld_[fi];
        }
    }

    volScalarField deltaRho(fvc::div(rhoPhi_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRho_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRho)
        );
    }
    deltaRho *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaRho += bi[fi]*deltaRho_[fi];
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    rho_.oldTime() = rhoOld;
    rho_ = rhoOld - dT*deltaRho;
    rho_.correctBoundaryConditions();

    thermo_->solve(stepi, ai, bi);
    phaseCompressibleSystem::solve(stepi, ai, bi);

    decode();
}


void Foam::singlePhaseCompressibleSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    phaseCompressibleSystem::setODEFields(nSteps, storeFields, storeDeltas);
    rhoOld_.setSize(nOld_);

    deltaRho_.setSize(nDelta_);
    thermo_->setODEFields(nSteps, oldIs_, nOld_, deltaIs_, nDelta_);
}

void Foam::singlePhaseCompressibleSystem::clearODEFields()
{
    phaseCompressibleSystem::clearODEFields();

    rhoOld_.clear();
    rhoOld_.setSize(nOld_);

    deltaRho_.clear();
    deltaRho_.setSize(nDelta_);
    thermo_->clearODEFields();
}


void Foam::singlePhaseCompressibleSystem::update()
{
    fluxScheme_->update
    (
        rho_,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::singlePhaseCompressibleSystem::ESource() const
{
    return thermo_->ESource();
}


void Foam::singlePhaseCompressibleSystem::decode()
{
    volScalarField rho(rho_);
    U_.ref() = rhoU_()/rho();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    e_.ref() = rhoE_()/rho() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();

    //- Limit internal energy it there is a negative temperature
    if(min(T_).value() < TLow_.value() && thermo_->limit())
    {
        if (debug)
        {
            WarningInFunction
                << "Lower limit of temperature reached, min(T) = "
                << min(T_).value()
                << ", limiting internal energy." << endl;
        }
        volScalarField limit(pos(T_ - TLow_));
        T_.max(TLow_);
        e_ = e_*limit + thermo_->E()*(1.0 - limit);
        rhoE_.ref() = rho_*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_->correct();
}


void Foam::singlePhaseCompressibleSystem::encode()
{
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volScalarField>
Foam::singlePhaseCompressibleSystem::speedOfSound() const
{
    return tmp<volScalarField>
    (
        new volScalarField("speedOfSound", thermo_->speedOfSound())
    );
}

// ************************************************************************* //
