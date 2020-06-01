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

#include "eThermoModel.H"

template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::volScalarFieldProperty
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


template<class BasicThermo, class ThermoType>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::cellSetProperty
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
Foam::eThermoModel<BasicThermo, ThermoType>::patchFieldProperty
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

template<class BasicThermo, class ThermoType>
Foam::eThermoModel<BasicThermo, ThermoType>::eThermoModel
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    BasicThermo
    (
        phaseName,
        p,
        rho,
        e,
        T,
        dict,
        master
    ),
    ThermoType(dict)
{}


template<class BasicThermo, class ThermoType>
Foam::eThermoModel<BasicThermo, ThermoType>::eThermoModel
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master
)
:
    BasicThermo
    (
        phaseName,
        mesh,
        dict,
        master
    ),
    ThermoType(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::eThermoModel<BasicThermo, ThermoType>::~eThermoModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::calcT() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("T", basicThermoModel::name_),
        dimTemperature,
        &ThermoType::TRhoE,
        this->T_,
        this->rho_,
        this->e_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::TRhoE
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
Foam::eThermoModel<BasicThermo, ThermoType>::TRhoEi
(
    const scalar& T,
    const scalar& e,
    const label celli
) const
{
    return ThermoType::TRhoE(T, this->rho_[celli], e);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::E() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("e", basicThermoModel::name_),
        dimEnergy/dimMass,
        &ThermoType::Es,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::e
(
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("e", basicThermoModel::name_),
        dimEnergy/dimMass,
        &ThermoType::Es,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::e
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
Foam::eThermoModel<BasicThermo, ThermoType>::e
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
Foam::eThermoModel<BasicThermo, ThermoType>::W() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("W", basicThermoModel::name_),
        dimMass/dimMoles,
        &ThermoType::W
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::W(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::W,
        patchi
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::eThermoModel<BasicThermo, ThermoType>::Wi(const label celli) const
{
    return ThermoType::W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::Gamma() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Gamma", basicThermoModel::name_),
        dimless,
        &ThermoType::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::Gamma(const label patchi) const
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
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::Cp() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cp", basicThermoModel::name_),
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cp,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::Cp(const label patchi) const
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
Foam::eThermoModel<BasicThermo, ThermoType>::Cp
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
Foam::eThermoModel<BasicThermo, ThermoType>::Cpi(const label celli) const
{
    return ThermoType::Cp(this->rho_[celli], this->e_[celli], this->T_[celli]);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::Cv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cv", basicThermoModel::name_),
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::Cv(const label patchi) const
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
Foam::eThermoModel<BasicThermo, ThermoType>::Cv
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
Foam::eThermoModel<BasicThermo, ThermoType>::Cvi(const label celli) const
{
    return ThermoType::Cv(this->rho_[celli], this->e_[celli], this->T_[celli]);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::CpByCv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("CpByCv", basicThermoModel::name_),
        dimless,
        &ThermoType::CpByCv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::eThermoModel<BasicThermo, ThermoType>::CpByCv(const label patchi) const
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
Foam::eThermoModel<BasicThermo, ThermoType>::CpByCv
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

// ************************************************************************* //
