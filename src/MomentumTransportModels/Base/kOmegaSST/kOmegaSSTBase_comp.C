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

#include "kOmegaSSTBase_comp.H"
#include "fluidThermo.H"
#include "fluidBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
void
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::boundOmega()
{
    omega_ = max(omega_, k_/(this->nutMaxCoeff_*this->nu()));
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        CDkOmegaMin_
        // dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*this->y()),
                scalar(500)*this->nu()/(sqr(this->y())*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(this->y()))
        ),
        arg1Max_
    );

    return tanh(pow4(arg1));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*this->y()),
            scalar(500)*this->nu()/(sqr(this->y())*omega_)
        ),
        arg2Max_
    );

    return tanh(sqr(arg2));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*this->nu()/(omega_*sqr(this->y())),
        arg3Max_
    );

    return 1 - tanh(pow4(arg3));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::
correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::
omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const HashTable<scalar>& defaults
)
:
    MomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    compressible::correction(this->coeffDict_, MODEL::K_OMEGA, *this),

    alphaK1_
    (
        lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            defaults,
            0.85
        )
    ),
    alphaK2_
    (
        lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            defaults,
            1.0
        )
    ),
    alphaOmega1_
    (
        lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            defaults,
            0.5
        )
    ),
    alphaOmega2_
    (
        lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            defaults,
            0.856
        )
    ),

    beta1_
    (
        lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            defaults,
            0.075
        )
    ),
    beta2_
    (
        lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            defaults,
            0.0828
        )
    ),
    betaStar_
    (
        lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            defaults,
            0.09
        )
    ),

    kappa_
    (
        lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            defaults,
            0.41
        )
    ),

    gamma1_
    (
        "gamma1",
        beta1_/betaStar_ - sqr(kappa_)*alphaOmega1_/sqrt(betaStar_)
    ),
    gamma2_
    (
        "gamma2",
        beta2_/betaStar_ - sqr(kappa_)*alphaOmega2_/sqrt(betaStar_)
    ),

    a1_
    (
        lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            defaults,
            0.31
        )
    ),
    b1_
    (
        lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            defaults,
            1.0
        )
    ),
    c1_
    (
        lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            defaults,
            20.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    arg1Max_
    (
        lookupOrAddToDict
        (
            "arg1Max",
            this->coeffDict_,
            defaults,
            10.0
        )
    ),
    arg2Max_
    (
        lookupOrAddToDict
        (
            "arg2Max",
            this->coeffDict_,
            defaults,
            10.0
        )
    ),
    arg3Max_
    (
        lookupOrAddToDict
        (
            "arg3Max",
            this->coeffDict_,
            defaults,
            10.0
        )
    ),
    CDkOmegaMin_
    (
        lookupOrAddToDict
        (
            "CDkOmegaMin",
            this->coeffDict_,
            defaults,
            dimless/sqr(dimTime),
            1.0e-10
        )
    ),

    k_
    (
        IOobject
        (
            this->groupName("k"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            this->groupName("omega"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    boundOmega();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
bool kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::read()
{
    if (MomentumTransportModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());

        kappa_.readIfPresent(this->coeffDict());

        gamma1_ =
            beta1_/betaStar_ - sqr(kappa_)*alphaOmega1_/sqrt(betaStar_);
        gamma2_ =
            beta2_/betaStar_ - sqr(kappa_)*alphaOmega2_/sqrt(betaStar_);

        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        arg1Max_.readIfPresent(this->coeffDict());
        arg2Max_.readIfPresent(this->coeffDict());
        arg3Max_.readIfPresent(this->coeffDict());
        CDkOmegaMin_.readIfPresent(this->coeffDict());

        compressible::correction::read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::correct()
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

    MomentumTransportModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)
       *(fvc::grad(this->k_) & fvc::grad(this->omega_))
       /this->omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNu,
                (this->c1_/this->a1_)*this->betaStar_*this->omega_()
               *max(this->a1_*this->omega_(), this->b1_*F23()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
          - fvm::SuSp
            (
                alpha()*rho()
               *(F1() - scalar(1))*CDkOmega()
               /this->omega_(),
                this->omega_
            )
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvModels.source(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        Foam::solve(omegaEqn);
        fvConstraints.constrain(this->omega_);
        this->boundOmega();
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
     ==
        alpha()*rho()*this->Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
      - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, F23), this->k_)
      + this->kSource()
      + fvModels.source(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    Foam::solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);
    this->boundOmega();

    this->correctNut(S2, F23);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
