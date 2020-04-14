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

#include "psiuCompressibleSystem.H"
#include "fiveEqnCompressibleTurbulenceModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiuCompressibleSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiuCompressibleSystem::psiuCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseCompressibleSystem(mesh, dict),
    thermo_(psiuReactionThermo::New(mesh)),
    e_(thermo_->he()),
    eu_(thermo_->heu()),
    rhoEu_
    (
        IOobject
        (
            "rhoEu",
            mesh.time().timeName(),
            mesh
        ),
        eu_*rho_,
        eu_.boundaryField().types()
    )
{
    rho_ = thermo_->rho();

    thermo_->validate("psiuCompressibleSystem", "ea");
    setModels(thermo_());

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiuCompressibleSystem::~psiuCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::psiuCompressibleSystem::solve
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
        rhoEuOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoEu_)
        );
    }

    volScalarField rhoOld(ai[stepi - 1]*rho_);
    volScalarField rhoEuOld(ai[stepi - 1]*rhoEu_);
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoOld += ai[fi]*rhoOld_[fi];
            rhoEuOld += ai[fi]*rhoEuOld_[fi];
        }
    }

    volScalarField deltaRho(fvc::div(rhoPhi_));
    volScalarField deltaRhoEu
    (
        fvc::div(fluxScheme_->energyFlux(rho_, U_, eu_, p_))
    );
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRho_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRho)
        );
        deltaRhoEu_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRhoEu)
        );
    }
    deltaRho *= bi[stepi - 1];
    deltaRhoEu *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaRho += bi[fi]*deltaRho_[fi];
            deltaRhoEu += bi[fi]*deltaRhoEu_[fi];
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    rho_ = rhoOld - dT*deltaRho;
    rho_.correctBoundaryConditions();
    rhoEu_ = rhoEuOld - dT*deltaRhoEu;

    phaseCompressibleSystem::solve(stepi, ai, bi);
    if (stepi == oldIs_.size() && turbulence_.valid())
    {
        Foam::solve
        (
            fvm::ddt(rho_, eu_) - fvc::ddt(rho_, eu_)
          - fvm::laplacian(turbulence_->alphaEff(), eu_)
        );
        rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));
        rhoEu_.boundaryFieldRef() =
            rho_.boundaryField()
           *(
                eu_.boundaryField()
              + 0.5*magSqr(U_.boundaryField())
            );
    }

//     if(temperatureFix)
    {
        scalar Tulow = 250.0;
        volScalarField dummyTu(thermo_->Tu());
        dummyTu.min(Tulow);

        eu_ = max(eu_, thermo_->he(p_, dummyTu));
        rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));
    }


    decode();
}


void Foam::psiuCompressibleSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    phaseCompressibleSystem::setODEFields(nSteps, storeFields, storeDeltas);
    rhoOld_.setSize(nOld_);
    rhoEuOld_.setSize(nOld_);

    deltaRho_.setSize(nDelta_);
    deltaRhoEu_.setSize(nDelta_);
}

void Foam::psiuCompressibleSystem::clearODEFields()
{
    phaseCompressibleSystem::clearODEFields();

    rhoOld_.clear();
    rhoOld_.setSize(nOld_);
    rhoEuOld_.clear();
    rhoEuOld_.setSize(nOld_);

    deltaRho_.clear();
    deltaRho_.setSize(nDelta_);
    deltaRhoEu_.clear();
    deltaRhoEu_.setSize(nDelta_);
}


void Foam::psiuCompressibleSystem::update()
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
Foam::psiuCompressibleSystem::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                rho_.mesh().time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho_.mesh(),
            dimensionedScalar("0", rhoE_.dimensions()/dimTime, 0.0)
        )
    );
}


void Foam::psiuCompressibleSystem::decode()
{
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/rho_);
    e_.ref() = E() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();

    volScalarField Eu(rhoEu_/rho_);
    eu_.ref() = Eu() - 0.5*magSqr(U_());
    eu_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );
    rhoEu_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            eu_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_->correct();
    p_.ref() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();
}


void Foam::psiuCompressibleSystem::encode()
{
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
    rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::speedOfSound() const
{
    return sqrt(thermo_->gamma()/thermo_->psi());
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::rhou() const
{
    return thermo_->rhou();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::Cv() const
{
    return thermo_->Cv();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::mu() const
{
    return thermo_->mu();
}


Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::mu(const label patchi) const
{
    return thermo_->mu(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::muu() const
{
    return thermo_->muu();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::mub() const
{
    return thermo_->mub();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::nu() const
{
    return thermo_->nu();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::nu(const label patchi) const
{
    return thermo_->nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::alpha() const
{
    return thermo_->alpha();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::alpha(const label patchi) const
{
    return thermo_->alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->alphaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::psiuCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->alphaEff(alphat, patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::alphahe() const
{
    return thermo_->alphahe();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::alphahe(const label patchi) const
{
    return thermo_->alphahe(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::kappa() const
{
    return thermo_->kappa();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::kappa(const label patchi) const
{
    return thermo_->kappa(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->kappaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::psiuCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->kappaEff(alphat, patchi);
}
// ************************************************************************* //
