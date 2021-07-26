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
    const fvMesh& mesh
)
:
    phaseCompressibleSystem(mesh),
    thermo_
    (
        fluidThermoModel::New
        (
            word::null,
            mesh,
            this->optionalSubDict("mixture"),
            true
        )
    )
{
    this->setModels();
    thermo_->initializeModels();
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::~singlePhaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singlePhaseCompressibleSystem::solve()
{
    volScalarField deltaRho("deltaRho", fvc::div(rhoPhi_));
    this->storeAndBlendDelta(deltaRho);

    dimensionedScalar dT = rho_.time().deltaT();
    this->storeAndBlendOld(rho_);
    rho_.storePrevIter();

    rho_ -= dT*deltaRho;
    rho_.correctBoundaryConditions();

    thermo_->solve();

    phaseCompressibleSystem::solve();
}


void Foam::singlePhaseCompressibleSystem::update()
{
    decode();
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
    thermo_->update();
}


Foam::tmp<Foam::volScalarField>
Foam::singlePhaseCompressibleSystem::ESource() const
{
    return thermo_->ESource();
}


void Foam::singlePhaseCompressibleSystem::decode()
{
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    e_.ref() = rhoE_()/rho_() - 0.5*magSqr(U_());

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
