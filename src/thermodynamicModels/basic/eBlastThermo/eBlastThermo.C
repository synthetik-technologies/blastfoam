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

#include "eBlastThermo.H"

template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::volScalarFieldProperty
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
            this->rho_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(this->rho_, celli)
    {
        psi[celli] = (this->*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->rho_.boundaryField()[patchi], facei)
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


template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::cellSetProperty
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


template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::patchFieldProperty
(
    Method psiMethod,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->rho_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->rho_.boundaryField()[patchi], facei)
    {
        psi[facei] =
            (this->*psiMethod)(args[facei] ...);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::eBlastThermo<BasicThermo, ThermoType>::eBlastThermo
(
    const word& phaseName,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
:
    BasicThermo
    (
        phaseName,
        rho,
        e,
        T,
        dict,
        masterName
    ),
    ThermoType(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::eBlastThermo<BasicThermo, ThermoType>::~eBlastThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("he", basicBlastThermo::name_),
        dimEnergy/dimMass,
        &ThermoType::Es,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::he
(
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &ThermoType::Es,
        faceCells,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::he
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Es,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("hs", basicBlastThermo::name_),
        dimEnergy/dimMass,
        &ThermoType::Hs,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::hs
(
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &ThermoType::Hs,
        faceCells,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Hs,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("ha", basicBlastThermo::name_),
        dimEnergy/dimMass,
        &ThermoType::Ha,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::ha
(
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &ThermoType::Ha,
        faceCells,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Ha,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::hc() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("hc", basicBlastThermo::name_),
        dimEnergy/dimMass,
        &ThermoType::Hf,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::flameT() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("flameT", basicBlastThermo::name_),
        dimTemperature,
        &ThermoType::flameT
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::THE() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("THE", basicBlastThermo::name_),
        dimTemperature,
        &ThermoType::TRhoE,
        T,
        this->rho_,
        this->e_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("THE", basicBlastThermo::name_),
        dimTemperature,
        &ThermoType::TRhoE,
        T,
        this->rho_,
        he
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &ThermoType::TRhoE,
        faceCells,
        T,
        this->rho_,
        he
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::TRhoE,
        patchi,
        T,
        this->rho_.boundaryField()[patchi],
        he
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::eBlastThermo<BasicThermo, ThermoType>::THEi
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return ThermoType::TRhoE(T, this->rho_[celli], he);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::Cp() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cp", basicBlastThermo::name_),
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cp,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Cp,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::eBlastThermo<BasicThermo, ThermoType>::Cpi
(
    const label celli
) const
{
    return
        ThermoType::Cp
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::Cv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cv", basicBlastThermo::name_),
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Cv,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::eBlastThermo<BasicThermo, ThermoType>::Cvi
(
    const label celli
) const
{
    return
        ThermoType::Cv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
}



template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::gamma() const
{
    return volScalarField::New("gamma", Cp()/Cv());
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp(T, patchi)/Cv(T, patchi);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::Cpv() const
{
    return volScalarField::New("Cpv", Cv());
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv(T, patchi);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::W() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("W", basicBlastThermo::name_),
        dimMass/dimMoles,
        &ThermoType::W
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::W(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::W,
        patchi
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::eBlastThermo<BasicThermo, ThermoType>::Wi(const label celli) const
{
    return ThermoType::W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::kappa() const
{
    return volScalarField::New
    (
        "kappa",
        Cp()*this->alpha_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::kappa
(
    const label patchi
) const
{
    return
        Cp(patchi)*gamma(this->T_[patchi], patchi)
       *this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::alphahe() const
{
    return volScalarField::New
    (
        "alphahe",
        gamma()*this->alpha_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::alphahe
(
    const label patchi
) const
{
    return
        gamma(this->T_[patchi], patchi)
       *this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "kappaEff",
        Cp()*(this->alpha_ + alphat)
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        Cp(gamma(this->T_[patchi], patchi), patchi)
       *(this->alpha_.boundaryField()[patchi] + alphat);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::alphaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "alphaEff",
        gamma()*(this->alpha_ + alphat)
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        gamma(gamma(this->T_[patchi], patchi), patchi)
       *(this->alpha_.boundaryField()[patchi] + alphat);
}


// ************************************************************************* //
