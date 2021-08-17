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
    )
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
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::initializeModels()
{
    activation_->initializeModels();
    afterburn_->initializeModels();

    Thermo::initializeModels();
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
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::speedOfSound() const
{
    tmp<volScalarField> cSqr
    (
        Thermo::blendedVolScalarFieldProperty
        (
            "cSqr",
            sqr(dimVelocity),
            &Thermo::thermoType1::cSqr,
            &Thermo::thermoType2::cSqr,
            Thermo::baseThermo::p(),
            this->rho_,
            this->e_,
            this->T_
        )
    );
    cSqr.ref().max(small);
    return volScalarField::New
    (
       this->phasePropertyName("speedOfSound"), sqrt(cSqr)
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidBlastThermo<Thermo>::speedOfSound(const label patchi) const
{
    return sqrt
    (
        max
        (
            Thermo::blendedPatchFieldProperty
            (
                &Thermo::thermoType1::cSqr,
                &Thermo::thermoType2::cSqr,
                patchi,
                Thermo::baseThermo::p().boundaryField()[patchi],
                this->rho_.boundaryField()[patchi],
                this->e_.boundaryField()[patchi],
                this->T_.boundaryField()[patchi]
            ),
            small
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::Gamma() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "Gamma",
        dimless,
        &Thermo::thermoType1::Gamma,
        &Thermo::thermoType2::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidBlastThermo<Thermo>::Gamma(const label patchi) const
{
    return Thermo::blendedPatchFieldProperty
    (
        &Thermo::thermoType1::Gamma,
        &Thermo::thermoType2::Gamma,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    const scalar& x = this->xi(celli);
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
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::pRhoT() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "p",
        dimPressure,
        &Thermo::thermoType1::p,
        &Thermo::thermoType2::p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellpRhoT(const label celli) const
{
    const scalar& x = this->xi(celli);
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
Foam::detonatingFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    const scalar& x = this->xi(celli);
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
    const scalar& x = this->xi(celli);
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
Foam::detonatingFluidBlastThermo<Thermo>::calcMu() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        &Thermo::thermoType1::mu,
        &Thermo::thermoType2::mu,
        this->rho_,
        this->e_,
        this->T_
    );
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
