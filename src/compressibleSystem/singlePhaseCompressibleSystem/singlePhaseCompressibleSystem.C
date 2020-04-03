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
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimVelocity), 0.0),
        fluidThermoModel::eBoundaryTypes(T_),
        fluidThermoModel::eBoundaryBaseTypes(T_)
    ),
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
    if (oldIs_[stepi - 1] != -1)
    {
        rhoOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rho_)
        );
    }

    volScalarField rhoOld(ai[stepi - 1]*rho_);
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

    //--- Hard limit, e
    if(min(e_).value() < 0)
    {
        WarningInFunction<< "Limiting e, min(e) = " << min(e_).value() << endl;
        e_.max(small);
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


Foam::tmp<Foam::volScalarField> Foam::singlePhaseCompressibleSystem::Cv() const
{
    return thermo_->Cv();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseCompressibleSystem::mu() const
{
    return thermo_->mu();
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseCompressibleSystem::mu(const label patchi) const
{
    return thermo_->mu(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseCompressibleSystem::nu() const
{
    return thermo_->nu();
}

Foam::tmp<Foam::scalarField>
Foam::singlePhaseCompressibleSystem::nu(const label patchi) const
{
    return thermo_->nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseCompressibleSystem::alpha() const
{
    return thermo_->alpha();
}

Foam::tmp<Foam::scalarField>
Foam::singlePhaseCompressibleSystem::alpha(const label patchi) const
{
    return thermo_->alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->alphaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::singlePhaseCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->alphaEff(alphat, patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseCompressibleSystem::alphahe() const
{
    return thermo_->alphahe();
}

Foam::tmp<Foam::scalarField>
Foam::singlePhaseCompressibleSystem::alphahe(const label patchi) const
{
    return thermo_->alphahe(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseCompressibleSystem::kappa() const
{
    return thermo_->kappa();
}

Foam::tmp<Foam::scalarField>
Foam::singlePhaseCompressibleSystem::kappa(const label patchi) const
{
    return
        thermo_->kappa(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->kappaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::singlePhaseCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        thermo_->kappaEff(alphat, patchi);
}
// ************************************************************************* //
