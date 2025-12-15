/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "kEpsilon_comp.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "fvc.H"
#include "fluxScheme.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilon_comp<BasicMomentumTransportModel>::boundEpsilon()
{
    tmp<volScalarField> tCmuk2(this->Cmu_*sqr(this->k_));
    this->epsilon_ = max
    (
        this->epsilon_,
        tCmuk2()/(this->nutMaxCoeff_*this->nu())
    );
    return tCmuk2;
}


template<class BasicMomentumTransportModel>
void kEpsilon_comp<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = boundEpsilon()/this->epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsilon_comp<BasicMomentumTransportModel>::kEpsilon_comp
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
    kEpsilon<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),
    compressible::correction(this->coeffDict_, MODEL::K_EPSILON, *this)
{
    bound(this->k_, this->kMin_);
    boundEpsilon();

    if (type == typeName || type == baseName())
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kEpsilon_comp<BasicMomentumTransportModel>::read()
{
    if (kEpsilon<BasicMomentumTransportModel>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kEpsilon_comp<BasicMomentumTransportModel>::correct()
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
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    this->epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, this->epsilon_)
      + fvm::div(alphaRhoPhi, this->epsilon_)
      - fvm::laplacian(alpha*rho*this->DepsilonEff(), this->epsilon_)
     ==
        this->C1_*alpha()*rho()*G*this->epsilon_()/this->k_()
      - fvm::SuSp
        (
            ((2.0/3.0)*this->C1_ - this->C3_)*alpha()*rho()*divU,
            this->epsilon_
        )
      - fvm::Sp
        (
            this->C2_*alpha()*rho()*this->epsilon_()/this->k_(),
            this->epsilon_
        )
      + this->epsilonSource()
      + fvModels.source(alpha, rho, this->epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(this->epsilon_.boundaryFieldRef());
    ::Foam::solve(epsEqn);
    fvConstraints.constrain(this->epsilon_);
    boundEpsilon();

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(), this->k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
      - fvm::Sp
        (
            alpha()*rho()*this->epsilon_()
           *(
                1.0/this->k_()
              + this->MtSqrByk()
            ),
            this->k_
        )
      + alpha()*rho()*this->pressureDialationSource(G)
      + this->kSource()
      + fvModels.source(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    ::Foam::solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
