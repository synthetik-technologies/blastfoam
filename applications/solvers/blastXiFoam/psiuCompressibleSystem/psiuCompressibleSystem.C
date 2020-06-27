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
    integrationSystem("phaseCompressibleSystem", mesh),
    thermo_(psiuReactionThermo::New(mesh)),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0.0),
        wordList(mesh.boundaryMesh().size(), "zeroGradient")
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_(thermo_->p()),
    T_(thermo_->T()),
    e_(thermo_->he()),
    eu_(thermo_->heu()),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho_*U_,
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity*sqr(dimVelocity), 0.0),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
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
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity*dimVelocity*dimArea, 0.0)
    ),
    rhoUPhi_
    (
        IOobject
        (
            "rhoUPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("0", dimDensity*sqr(dimVelocity)*dimArea, Zero)
    ),
    rhoEPhi_
    (
        IOobject
        (
            "rhoEPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity*pow3(dimVelocity)*dimArea, 0.0)
    ),
    fluxScheme_(fluxScheme::New(mesh)),
    g_(mesh.lookupObject<uniformDimensionedVectorField>("g"))
{
    rho_ = thermo_->rho();

    thermo_->validate("psiuCompressibleSystem", "ea");

    if (min(thermo_->mu()).value() > small)
    {
        turbulence_.set
        (
            compressible::turbulenceModel::New
            (
                rho_,
                U_,
                rhoPhi_,
                thermo_()
            ).ptr()
        );
    }

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
        rhoUOld_.set
        (
            oldIs_[stepi - 1],
            new volVectorField(rhoU_)
        );
        rhoEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoE_)
        );
        rhoEuOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoEu_)
        );
    }

    volScalarField rhoOld(ai[stepi - 1]*rho_);
    volVectorField rhoUOld(ai[stepi - 1]*rhoU_);
    volScalarField rhoEOld(ai[stepi - 1]*rhoE_);
    volScalarField rhoEuOld(ai[stepi - 1]*rhoEu_);
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoOld += ai[fi]*rhoOld_[fi];
            rhoUOld += ai[fi]*rhoUOld_[fi];
            rhoEOld += ai[fi]*rhoEOld_[fi];
            rhoEuOld += ai[fi]*rhoEuOld_[fi];
        }
    }

    volScalarField deltaRho(fvc::div(rhoPhi_));
    volVectorField deltaRhoU(fvc::div(rhoUPhi_) - g_*rho_);
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - (rhoU_ & g_)
    );
    volScalarField deltaRhoEu
    (
        fvc::div(fluxScheme_->energyFlux(rho_, U_, eu_, p_))
      - (rhoU_ & g_)
    );
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRho_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRho)
        );
        deltaRhoU_.set
        (
            deltaIs_[stepi - 1],
            new volVectorField(deltaRhoU)
        );
        deltaRhoE_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRhoE)
        );
        deltaRhoEu_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRhoEu)
        );
    }
    deltaRho *= bi[stepi - 1];
    deltaRhoU *= bi[stepi - 1];
    deltaRhoE *= bi[stepi - 1];
    deltaRhoEu *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaRho += bi[fi]*deltaRho_[fi];
            deltaRhoU += bi[fi]*deltaRhoU_[fi];
            deltaRhoE += bi[fi]*deltaRhoE_[fi];
            deltaRhoEu += bi[fi]*deltaRhoEu_[fi];
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    rho_ = rhoOld - dT*deltaRho;
    rho_.correctBoundaryConditions();

    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);
    rhoU_ = cmptMultiply(rhoUOld - dT*deltaRhoU, solutionDs);
    rhoE_ = rhoEOld - dT*deltaRhoE;
    rhoEu_ = rhoEuOld - dT*deltaRhoEu;

    if (stepi == oldIs_.size() && (turbulence_.valid()))
    {
        U_ = rhoU_/rho_;
        U_.correctBoundaryConditions();

        volScalarField muEff("muEff", turbulence_->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U_))));

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
         ==
            fvm::laplacian(muEff, U_)
          + fvc::div(tauMC)
        );

        UEqn.solve();

        rhoU_ = rho_*U_;

        e_ = rhoE_/rho_ - 0.5*magSqr(U_);
        e_.correctBoundaryConditions();

        eu_ = rhoEu_/rho_ - 0.5*magSqr(U_);
        eu_.correctBoundaryConditions();

        Foam::solve
        (
            fvm::ddt(rho_, e_) - fvc::ddt(rho_, e_)
        - fvm::laplacian(turbulence_->alphaEff(), e_)
        );
        rhoE_ = rho_*(e_ + 0.5*magSqr(U_)); // Includes change to total energy from viscous term in momentum equation

        Foam::solve
        (
            fvm::ddt(rho_, eu_) - fvc::ddt(rho_, eu_)
        - fvm::laplacian(turbulence_->alphaEff(), eu_)
        );
        rhoE_ = rho_*(eu_ + 0.5*magSqr(U_)); // Includes change to total energy from viscous term in momentum equation

        turbulence_->correct();
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
    oldIs_.resize(nSteps);
    deltaIs_.resize(nSteps);
    label fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            oldIs_[i] = fi;
            fi++;
        }
        else
        {
            oldIs_[i] = -1;
        }
    }
    nOld_ = fi;
    rhoOld_.resize(nOld_);
    rhoUOld_.resize(nOld_);
    rhoEOld_.resize(nOld_);
    rhoEuOld_.resize(nOld_);

    fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaIs_[i] = fi;
            fi++;
        }
        else
        {
            deltaIs_[i] = -1;
        }
    }
    nDelta_ = fi;
    deltaRho_.resize(nDelta_);
    deltaRhoU_.resize(nDelta_);
    deltaRhoE_.resize(nDelta_);
    deltaRhoEu_.resize(nDelta_);
}

void Foam::psiuCompressibleSystem::clearODEFields()
{
    fluxScheme_->clear();
    rhoOld_.clear();
    rhoUOld_.clear();
    rhoEOld_.clear();
    rhoEuOld_.clear();
    rhoOld_.resize(nOld_);
    rhoUOld_.resize(nOld_);
    rhoEOld_.resize(nOld_);
    rhoEuOld_.resize(nOld_);

    deltaRho_.clear();
    deltaRhoU_.clear();
    deltaRhoE_.clear();
    deltaRhoEu_.clear();
    deltaRho_.resize(nDelta_);
    deltaRhoU_.resize(nDelta_);
    deltaRhoE_.resize(nDelta_);
    deltaRhoEu_.resize(nDelta_);
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


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::rho() const
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::rhou() const
{
    return thermo_->rhou();
}


Foam::tmp<Foam::volVectorField> Foam::psiuCompressibleSystem::U() const
{
    return U_;
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
