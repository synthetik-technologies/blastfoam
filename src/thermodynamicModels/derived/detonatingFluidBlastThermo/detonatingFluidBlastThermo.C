/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a detonating material
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

#include "detonatingFluidBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::HE
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Es(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Es(rho, e, T);
    }
    return thermo2_.Es(rho, e, T)*xi_ + thermo1_.Es(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::TRhoE
(
    const scalar T,
    const scalar rho,
    const scalar e
) const
{
    if (xi_ < small)
    {
        return thermo1_.TRhoE(T, rho, e);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.TRhoE(T, rho, e);
    }

    return
        thermo2_.TRhoE(T, rho, e)*xi_
      + thermo1_.TRhoE(T, rho, e)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Cp(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Cp(rho, e, T);
    }
    return thermo2_.Cp(rho, e, T)*xi_ + thermo1_.Cp(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Cv(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Cv(rho, e, T);
    }
    return thermo2_.Cv(rho, e, T)*xi_ + thermo1_.Cv(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.kappa(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.kappa(rho, e, T);
    }
    return thermo2_.kappa(rho, e, T)*xi_ + thermo1_.kappa(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::pRhoT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.p(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.p(rho, e, T);
    }
    return thermo2_.p(rho, e, T)*xi_ + thermo1_.p(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Gamma(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Gamma(rho, e, T);
    }
    return thermo2_.Gamma(rho, e, T)*xi_ + thermo1_.Gamma(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::mu
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.mu(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.mu(rho, e, T);
    }
    return thermo2_.mu(rho, e, T)*xi_ + thermo1_.mu(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::speedOfSound
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return sqrt(max(cSqr(p, rho, e, T), small));
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastFluidMixture<BasicMixture, ThermoType1, ThermoType2>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.cSqr(p, rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.cSqr(p, rho, e, T);
    }
    return
        thermo2_.cSqr(p, rho, e, T)*xi_
      + thermo1_.cSqr(p, rho, e, T)*(1.0 - xi_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingFluidBlastThermo<Thermo>::detonatingFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    Thermo
    (
        mesh,
        dict,
        dict.subDict("reactants"),
        dict.subDict("products"),
        phaseName,
        masterName
    ),
    activation_
    (
        activationModel::New
        (
            mesh,
            dict,
            phaseName
        )
    ),
    afterburn_
    (
        afterburnModel::New
        (
            mesh,
            dict,
            phaseName
        )
    ),
    mixture_(*this, *this, activation_())
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() <= 0
     || (
            Thermo::thermoType1::solid()
         && dict.lookupOrDefault<Switch>("calculateDensity", false)
         && !mesh.time().restart()
        )
    )
    {
        updateRho(Thermo::baseThermo::p());
    }
    this->initializeFields();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::initializeModels()
{
    activation_->initializeModels();
    afterburn_->initializeModels();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingFluidBlastThermo<Thermo>::~detonatingFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::update()
{
    activation_->update();
    afterburn_->update();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::solve()
{
    activation_->solve();
    afterburn_->solve();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::postUpdate()
{
    activation_->postUpdate();
    afterburn_->postUpdate();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::updateRho(const volScalarField& p)
{
    volScalarField rhoNew
    (
        Thermo::blendedVolScalarFieldProperty
        (
            "rho",
            dimDensity,
            &Thermo::thermoType1::rhoPT,
            &Thermo::thermoType2::rhoPT,
            p,
            this->T_
        )
    );
    this->rho_ = rhoNew;
    forAll(this->rho_.boundaryField(), patchi)
    {
        forAll(this->rho_.boundaryField()[patchi], facei)
        {
            this->rho_.boundaryFieldRef()[patchi][facei] =
                rhoNew.boundaryField()[patchi][facei];
        }
    }
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellpRhoT(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::p(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::p(rho, e, T);
    }

    return
        Thermo::thermoType2::p(rho, e, T)*x
      + Thermo::thermoType1::p(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::Gamma(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::Gamma(rho, e, T);
    }

    return
        Thermo::thermoType2::Gamma(rho, e, T)*x
      + Thermo::thermoType1::Gamma(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::dpdRho(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::dpdRho(rho, e, T);
    }

    return
        Thermo::thermoType2::dpdRho(rho, e, T)*x
      + Thermo::thermoType1::dpdRho(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::celldpde(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::dpde(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::dpde(rho, e, T);
    }

    return
        Thermo::thermoType2::dpde(rho, e, T)*x
      + Thermo::thermoType1::dpde(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::calce(const volScalarField& p) const
{
    tmp<volScalarField> eInit
    (
        Thermo::volScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &Thermo::thermoType1::initializeEnergy,
            p,
            this->rho_,
            this->e_,
            this->T_
        )
    );

    //- Add detonation energy to initially reacted material
    if (this->rho_.time().timeIndex() == 0)
    {
        eInit.ref() += activation_->e0()*activation_->lambda();
    }

    return eInit;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (activation_->ESource() + afterburn_->ESource())*this->rho_
        )
    );
}



// ************************************************************************* //
