/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "kOmega_comp.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmega_comp<BasicMomentumTransportModel>::kOmega_comp
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
    kOmega<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),
    ::Foam::compressible::correction(this->coeffDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmega_comp<BasicMomentumTransportModel>::read()
{
    if (kOmega<BasicMomentumTransportModel>::read())
    {
        return ::Foam::compressible::correction::read(this->coeffDict_);
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kOmega_comp<BasicMomentumTransportModel>::correct()
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

    tmp<volScalarField::Internal> beta, betaStar;
    ::Foam::compressible::correction::correct
    (
        this->k_,
        this->beta_,
        this->betaStar_,
        beta,
        betaStar
    );

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

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

    this->correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
