/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun     : Derived from compressibilityCorrection
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "kOmega2006_comp.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kOmega2006_comp<BasicMomentumTransportModel>::correctNut
(
    const volTensorField& gradU
)
{
    this->nut_ =
        this->k_
       /max(this->omega_, this->Clim_*sqrt(2/this->betaStar_)*mag(dev(symm(gradU))));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> kOmega2006_comp<BasicMomentumTransportModel>::beta
(
    const volTensorField& gradU
) const
{
    const volSymmTensorField::Internal S(symm(gradU()));
    const volSymmTensorField::Internal Shat(S - 0.5*tr(S)*I);
    const volTensorField::Internal Omega(skew(gradU.v()));

    const volScalarField::Internal ChiOmega
    (
        typedName("ChiOmega"),
        mag((Omega & Omega) && Shat)/pow3(this->betaStar_*this->omega_.v())
    );

    const volScalarField::Internal fBeta
    (
        typedName("fBeta"),
        (1 + 85*ChiOmega)/(1 + 100*ChiOmega)
    );

    return this->beta0_*fBeta;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmega2006_comp<BasicMomentumTransportModel>::CDkOmega() const
{
    return max
    (
        this->sigmaDo_
       *(fvc::grad(this->k_)().v() & fvc::grad(this->omega_)().v())
       /this->omega_(),
        dimensionedScalar(dimless/sqr(dimTime), 0)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmega2006_comp<BasicMomentumTransportModel>::kOmega2006_comp
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    kOmega2006<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),
    ::Foam::compressible::correction(this->coeffDict_, SARKAR)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmega2006_comp<BasicMomentumTransportModel>::read()
{
    if (kOmega2006<BasicMomentumTransportModel>::read())
    {
        return ::Foam::compressible::correction::read(this->coeffDict_);
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kOmega2006_comp<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        typedName("divU"),
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    const volTensorField gradU(fvc::grad(U));

    tmp<volScalarField::Internal> beta(this->beta(gradU));
    tmp<volScalarField::Internal> betaStar;
    Foam::compressible::correction::correct
    (
        this->k_,
        this->betaStar_,
        beta,
        betaStar
    );

    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(gradU.v())) && gradU.v())
    );

    ::Foam::compressible::correction::limitG
    (
        G,
        this->k_,
        this->omega_,
        betaStar
    );

    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, this->omega_)
      + fvm::div(alphaRhoPhi, this->omega_)
      - fvm::laplacian(alpha*rho*this->DomegaEff(), this->omega_)
     ==
        this->gamma_*alpha()*rho()*G*this->omega_()/this->k_()
      - fvm::SuSp(((2.0/3.0)*this->gamma_)*alpha()*rho()*divU, this->omega_)
      - fvm::Sp(beta*alpha()*rho()*this->omega_(), this->omega_)
      + alpha()*rho()*this->CDkOmega()
      + this->omegaSource()
      + fvModels.source(alpha, rho, this->omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvConstraints.constrain(this->omega_);
    this->boundOmega();


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(), this->k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
      - fvm::Sp(betaStar*alpha()*rho()*this->omega_(), this->k_)
      + this->kSource()
      + fvModels.source(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);
    this->boundOmega();

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
