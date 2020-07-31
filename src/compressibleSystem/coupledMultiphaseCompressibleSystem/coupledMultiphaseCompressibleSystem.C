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

#include "coupledMultiphaseCompressibleSystem.H"
#include "blastCompressibleTurbulenceModel.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledMultiphaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        phaseCompressibleSystem,
        coupledMultiphaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMultiphaseCompressibleSystem::coupledMultiphaseCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    multiphaseCompressibleSystem(mesh, dict),
    volumeFraction_
    (
        IOobject
        (
            "continuousVolumeFraction",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        0.0
    ),
    alphaRho_
    (
        IOobject
        (
            "alphaRho",
            mesh.time().timeName(),
            mesh
        ),
        volumeFraction_*rho_
    )
{
    forAll(alphas_, phasei)
    {
        volumeFraction_ += alphas_[phasei];
    }
    this->thermo_.setTotalVolumeFractionPtr(volumeFraction_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledMultiphaseCompressibleSystem::~coupledMultiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMultiphaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    PtrList<volScalarField> alphasOld(alphas_.size());
    PtrList<volScalarField> alphaRhosOld(alphas_.size());
    forAll(alphas_, phasei)
    {
        alphasOld.set
        (
            phasei, new volScalarField(alphas_[phasei])
        );
        this->storeOld(stepi, alphasOld[phasei], alphasOld_[phasei]);
        this->blendOld(stepi, alphasOld[phasei], alphasOld_[phasei], ai);

        alphaRhosOld.set
        (
            phasei, new volScalarField(alphaRhos_[phasei])
        );
        this->storeOld(stepi, alphaRhosOld[phasei], alphaRhosOld_[phasei]);
        this->blendOld(stepi, alphaRhosOld[phasei], alphaRhosOld_[phasei], ai);
        alphaRhos_[phasei].oldTime() = alphaRhosOld[phasei];
    }

    rho_.oldTime() = alphaRhosOld[0];
    for (label i = 1; i < alphaRhosOld.size(); i++)
    {
        rho_.oldTime() += alphaRhosOld[i];
    }
    rho_.oldTime() /= max(volumeFraction_, 1e-10);

    volVectorField rhoUOld(rhoU_);
    this->storeOld(stepi, rhoUOld, rhoUOld_);
    this->blendOld(stepi, rhoUOld, rhoUOld_, ai);

    volScalarField rhoEOld(rhoE_);
    this->storeOld(stepi, rhoEOld, rhoEOld_);
    this->blendOld(stepi, rhoEOld, rhoEOld_, ai);

    thermo_.solve(stepi, ai, bi);

    surfaceScalarField alphaf(fluxScheme_->interpolate(volumeFraction_, "alpha"));
    PtrList<volScalarField> deltaAlphas(alphas_.size());
    PtrList<volScalarField> deltaAlphaRhos(alphas_.size());
    forAll(alphas_, phasei)
    {
        deltaAlphas.set
        (
            phasei,
            new volScalarField
            (
                fvc::div(alphaPhis_[phasei])
              - alphas_[phasei]*fvc::div(phi_)
            )
        );
        this->storeDelta(stepi, deltaAlphas[phasei], deltaAlphas_[phasei]);
        this->blendDelta(stepi, deltaAlphas[phasei], deltaAlphas_[phasei], bi);

        deltaAlphaRhos.set
        (
            phasei, new volScalarField(fvc::div(alphaRhoPhis_[phasei]))
        );
        this->storeDelta(stepi, deltaAlphaRhos[phasei], deltaAlphaRhos_[phasei]);
        this->blendDelta
        (
            stepi,
            deltaAlphaRhos[phasei],
            deltaAlphaRhos_[phasei],
            bi
        );
    }

    volVectorField deltaRhoU
    (
        fvc::div(rhoUPhi_*alphaf)
      - p_*fvc::grad(alphaf)
      - g_*alphaRho_
    );
    this->storeDelta(stepi, deltaRhoU, deltaRhoU_);
    this->blendDelta(stepi, deltaRhoU, deltaRhoU_, bi);

    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_*alphaf)
      - ESource()
      - (rhoU_ & g_)
    );
    this->storeDelta(stepi, deltaRhoE, deltaRhoE_);
    this->blendDelta(stepi, deltaRhoE, deltaRhoE_, bi);

    dimensionedScalar dT = rho_.time().deltaT();
    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);

    rhoU_ = cmptMultiply(rhoUOld - dT*deltaRhoU, solutionDs);
    rhoE_ = rhoEOld - dT*deltaRhoE;
    if (radiation_->type() != "none")
    {
        NotImplemented;
//         calcAlphaAndRho();
//         e() = rhoE_/alphaRho_ - 0.5*magSqr(U_);
//         e().correctBoundaryConditions();
//         rhoE_ = radiation_->calcRhoE(f*dT, rhoE_, alphaRho_, e(), Cv());
    }

    forAll(alphas_, phasei)
    {
        alphas_[phasei] = alphasOld[phasei] - dT*deltaAlphas[phasei];
        alphas_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei].oldTime() = alphaRhosOld[phasei];
        alphaRhos_[phasei] = alphaRhosOld[phasei] - dT*deltaAlphaRhos[phasei];
    }
}

void Foam::coupledMultiphaseCompressibleSystem::calcAlphaAndRho()
{
    volumeFraction_ = min(1.0, max(0.0, 1.0 - *alphadPtr_));
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            rho_.time().timeName(),
            rho_.mesh()
        ),
        rho_.mesh(),
        0.0
    );
    for (label phasei = 0; phasei < alphas_.size() - 1; phasei++)
    {
        alphas_[phasei].max(0);
        alphas_[phasei].min(1);
        alphas_[phasei].correctBoundaryConditions();
        sumAlpha += alphas_[phasei];

        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());
        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }
    label phasei = alphas_.size() - 1;
    alphas_[phasei] = volumeFraction_ - sumAlpha;
    alphas_[phasei].max(0.0);
    alphas_[phasei].min(1.0);

    alphaRhos_[phasei].max(0.0);
    rhos_[phasei] = alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());
    rhos_[phasei].correctBoundaryConditions();
    alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
    alphaRho_ += alphaRhos_[phasei];
    rho_ = alphaRho_/max(volumeFraction_, 1e-10);
}


void Foam::coupledMultiphaseCompressibleSystem::decode()
{
    calcAlphaAndRho();

    volScalarField alphaRhos(alphaRho_);
    alphaRhos.max(1e-10);
    U_.ref() = rhoU_()/alphaRhos();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = alphaRho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/alphaRhos);
    e_.ref() = E() - 0.5*magSqr(U_());

    //- Limit internal energy it there is a negative temperature
    if(min(T_).value() < TLow_.value() && thermo_.limit())
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
        e_ = e_*limit + thermo_.E()*(1.0 - limit);
        rhoE_.ref() = rho_*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        alphaRho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_.correct();
    T_ /= Foam::max(volumeFraction_, 1e-10);
}


void Foam::coupledMultiphaseCompressibleSystem::encode()
{
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }
    rho_ = alphaRho_/max(volumeFraction_, 1e-10);
    rhoU_ = alphaRho_*U_;
    rhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}



Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::Cv() const
{
    return thermo_.Cv();
}


Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::mu() const
{
    return thermo_.mu();
}


Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::mu(const label patchi) const
{
    return thermo_.mu(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::nu() const
{
    return thermo_.nu();
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::nu(const label patchi) const
{
    return thermo_.nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::coupledMultiphaseCompressibleSystem::alpha() const
{
    return thermo_.alpha();
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::alpha(const label patchi) const
{
    return thermo_.alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_.alphaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::coupledMultiphaseCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_.alphaEff(alphat, patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::coupledMultiphaseCompressibleSystem::alphahe() const
{
    return thermo_.alphahe();
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::alphahe(const label patchi) const
{
    return thermo_.alphahe(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::kappa() const
{
    return thermo_.kappa();
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::kappa(const label patchi) const
{
    return thermo_.kappa(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_.kappaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::coupledMultiphaseCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_.kappaEff(alphat, patchi);
}

// ************************************************************************* //
