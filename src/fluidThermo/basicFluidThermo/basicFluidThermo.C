/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "basicFluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Thermo>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::volScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            this->p_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(this->p_, celli)
    {
        psi[celli] = (this->*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                (this->*psiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                );
        }
    }

    return tPsi;
}


template<class Thermo>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField> Foam::basicFluidThermo<Thermo>::cellSetProperty
(
    Method psiMethod,
    const labelList& cells,
    const Args& ... args
) const
{
    // Note: Args are fields for the set, not for the mesh as a whole. The
    // cells list is only used to get the mixture.

    tmp<scalarField> tPsi(new scalarField(cells.size()));
    scalarField& psi = tPsi.ref();

    forAll(cells, celli)
    {
        psi[celli] =
           (this->*psiMethod)(args[celli] ...);
    }

    return tPsi;
}


template<class Thermo>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField> Foam::basicFluidThermo<Thermo>::patchFieldProperty
(
    Method psiMethod,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->p_.boundaryField()[patchi], facei)
    {
        psi[facei] =
            (this->*psiMethod)(args[facei] ...);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidThermo<Thermo>::basicFluidThermo
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    fluidThermoModel
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master
    ),
    Thermo(dict)
{
    if (dict.lookupOrDefault<Switch>("calculateDensity", false))
    {
        volScalarField rhoInit
        (
            volScalarFieldProperty
            (
                IOobject::groupName("rhoInit", name_),
                dimDensity,
                &Thermo::initializeRho,
                p_,
                rho_,
                e_,
                T_
            )
        );
        rho_ = rhoInit;
        forAll(rho_.boundaryField(), patchi)
        {
            forAll(rho_.boundaryField()[patchi], facei)
            {
                rho_.boundaryFieldRef()[patchi][facei] =
                    rhoInit.boundaryField()[patchi][facei];
            }
        }
    }
    if (master)
    {
        volScalarField e(calce());
        e_ = e;
        forAll(e_.boundaryField(), patchi)
        {
            forAll(e_.boundaryField()[patchi], facei)
            {
                e_.boundaryFieldRef()[patchi][facei] =
                    e.boundaryField()[patchi][facei];
            }
        }
        correct();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidThermo<Thermo>::~basicFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::basicFluidThermo<Thermo>::correct()
{
    if (master_)
    {
        T_ = calcT();
        p_ = calcP();
    }

    mu_ = volScalarFieldProperty
    (
        IOobject::groupName("mu", name_),
        dimDynamicViscosity,
        &Thermo::mu,
        rho_,
        e_,
        T_
    );

    alpha_ = volScalarFieldProperty
    (
        IOobject::groupName("alphah", name_),
        dimensionSet(1, -1, -1, 0, 0),
        &Thermo::alphah,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::speedOfSound() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("speedOfSound", name_),
        dimVelocity,
        &Thermo::speedOfSound,
        p_,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::speedOfSound(const label patchi) const
{
    return patchFieldProperty
    (
        &Thermo::speedOfSound,
        patchi,
        p_.boundaryField()[patchi],
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::calcT() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("T", name_),
        dimTemperature,
        &Thermo::TRhoE,
        T_,
        rho_,
        e_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::calcP() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("p", name_),
        dimPressure,
        &Thermo::p,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::calce() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("e", name_),
        dimEnergy/dimMass,
        &Thermo::initializeEnergy,
        p_,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::E() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("e", name_),
        dimEnergy/dimMass,
        &Thermo::Ea,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::E
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &Thermo::Ea,
        patchi,
        rho,
        e,
        T
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::E
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &Thermo::Ea,
        faceCells,
        rho,
        e,
        T
    );
}

template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::Gamma() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Gamma", name_),
        dimless,
        &Thermo::Gamma,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::Gamma(const label patchi) const
{
    return patchFieldProperty
    (
        &Thermo::Gamma,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("ESource", name_),
                p_.mesh().time().timeName(),
                p_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p_.mesh(),
            dimensionedScalar(sqr(dimVelocity)*dimDensity/dimTime, 0.0)
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::W() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("W", name_),
        dimMass/dimMoles,
        &Thermo::W
    );
}


template<class Thermo>
Foam::scalar Foam::basicFluidThermo<Thermo>::Wi(const label celli) const
{
    return Thermo::W();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::Cp() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cp", name_),
        dimEnergy/dimMass/dimTemperature,
        &Thermo::Cp,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::Cp(const label patchi) const
{
    return patchFieldProperty
    (
        &Thermo::Cp,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::scalar Foam::basicFluidThermo<Thermo>::Cpi(const label celli) const
{
    return Thermo::Cp(rho_[celli], e_[celli], T_[celli]);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::Cv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cv", name_),
        dimEnergy/dimMass/dimTemperature,
        &Thermo::Cv,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::Cv(const label patchi) const
{
    return patchFieldProperty
    (
        &Thermo::Cv,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &Thermo::Cv,
        patchi,
        rho,
        e,
        T
    );
}


template<class Thermo>
Foam::scalar Foam::basicFluidThermo<Thermo>::Cvi(const label celli) const
{
    return Thermo::Cv(rho_[celli], e_[celli], T_[celli]);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidThermo<Thermo>::CpByCv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("CpByCv", name_),
        dimless,
        &Thermo::CpByCv,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::CpByCv(const label patchi) const
{
    return patchFieldProperty
    (
        &Thermo::CpByCv,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}

// ************************************************************************* //
