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

#include "mixtureThermoModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
const Foam::PtrList<ThermoType>&
Foam::mixtureThermoModel<BasicThermo, ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(this->species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(this->species_[i]))
        );
    }

    return speciesData_;
}

template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::volScalarFieldProperty
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
        psi[celli] = (this->mixture_[celli].*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                (this->mixture_.boundary(patchi, facei).*psiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                );
        }
    }

    return tPsi;
}


template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::volScalarFieldSpecieProperty
(
    const word& s,
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

    const ThermoType& thermo = speciesData_[this->species_[s]];

    volScalarField& psi = tPsi.ref();

    forAll(this->p_, celli)
    {
        psi[celli] = (thermo.*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                (thermo.*psiMethod)
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
Foam::mixtureThermoModel<BasicThermo, ThermoType>::patchFieldProperty
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
    forAll(psi, facei)
    {
        psi[facei] =
            (this->mixture_.boundary(patchi, facei).*psiMethod)
            (
                args[facei] ...
            );
    }

    return tPsi;
}


template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::cellSetProperty
(
    Method psiMethod,
    const labelList& cells,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi(new scalarField(cells.size()));

    scalarField& psi = tPsi.ref();
    forAll(cells, celli)
    {
        psi[celli] =
            (this->mixture_[cells[celli]].*psiMethod)(args[celli] ...);
    }

    return tPsi;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::mixtureThermoModel
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    BasicThermo
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    ),
    speciesData_(this->species_.size()),
    mixture_(p.mesh(), this->Ys_, constructSpeciesData(dict), name)
{}


template<class BasicThermo, class ThermoType>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::mixtureThermoModel
(
    const HashPtrTable<ThermoType, word, string::hash>& thermoData,
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    BasicThermo
    (
        thermoData.toc(),
        name,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    ),
    speciesData_(this->species_.size()),
    mixture_(p.mesh(), this->Ys_, constructSpeciesData(dict), name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::~mixtureThermoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::calcT() const
{
    return volScalarFieldProperty
    (
        "T",
        dimTemperature,
        &ThermoType::TRhoE,
        this->T_,
        this->rho_,
        this->e_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::TRhoE
(
    const scalarField& T,
    const scalarField& e,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::TRhoE,
        patchi,
        T,
        this->rho_.boundaryField()[patchi],
        e
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::TRhoEi
(
    const scalar& T,
    const scalar& e,
    const label celli
) const
{
    return mixture_[celli].TRhoE(T, this->rho_[celli], e);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::E() const
{
    return volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &ThermoType::Es,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::e
(
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &ThermoType::Es,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Es,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &ThermoType::Es,
        faceCells,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::W() const
{
    return volScalarFieldProperty
    (
        "W",
        dimMass/dimMoles,
        &ThermoType::W
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::W(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::W,
        patchi
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Wi(const label celli) const
{
    return mixture_[celli].W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Gamma() const
{
    return volScalarFieldProperty
    (
        "Gamma",
        dimless,
        &ThermoType::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Gamma(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::Gamma,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Gammai(const label celli) const
{
    return mixture_[celli].Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cp() const
{
    return volScalarFieldProperty
    (
        "Cp",
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cp,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cp(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::Cp,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Cp,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cpi(const label celli) const
{
    return mixture_[celli].Cp
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cv() const
{
    return volScalarFieldProperty
    (
        "Cv",
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cv(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::Cv,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::Cv,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Cvi(const label celli) const
{
    return mixture_[celli].Cv
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::CpByCv() const
{
    return volScalarFieldProperty
    (
        "CpByCv",
        dimless,
        &ThermoType::CpByCv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::CpByCv(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::CpByCv,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::CpByCv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &ThermoType::CpByCv,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::Hf() const
{
    return volScalarFieldProperty
    (
        "Hf",
        dimEnergy/dimMass,
        &ThermoType::Hf
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::flameT() const
{
    return volScalarFieldProperty
    (
        "flameT",
        dimTemperature,
        &ThermoType::flameT
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieW(const word& s) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "W",
        dimMass/dimMoles,
        &ThermoType::W
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieWi
(
    const word& s,
    const label celli
) const
{
    return speciesData_[this->species_[s]].W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieGamma
(
    const word& s
) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "Gamma",
        dimless,
        &ThermoType::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieGammai
(
    const word& s,
    const label celli
) const
{
    return speciesData_[this->species_[s]].Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieCp(const word& s) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "Cp",
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cp,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieCpi
(
    const word& s,
    const label celli
) const
{
    return speciesData_[this->species_[s]].Cp
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieCv(const word& s) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "Cv",
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieCvi
(
    const word& s,
    const label celli
) const
{
    return speciesData_[this->species_[s]].Cv
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}

template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieCpByCv
(
    const word& s
) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "CpByCv",
        dimless,
        &ThermoType::CpByCv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieHf(const word& s) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "Hf",
        dimEnergy/dimMass,
        &ThermoType::Hf
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureThermoModel<BasicThermo, ThermoType>::specieFlameT
(
    const word& s
) const
{
    return volScalarFieldSpecieProperty
    (
        s,
        "flameT",
        dimTemperature,
        &ThermoType::flameT
    );
}

// ************************************************************************* //
