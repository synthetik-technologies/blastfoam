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
            IOobject::groupName(psiName, this->phaseName_),
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
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    BasicThermo
    (
        mesh,
        dict,
        phaseName,
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
Foam::word Foam::eBlastThermo<BasicThermo, ThermoType>::thermoName() const
{
    return ThermoType::typeName();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::calce
(
    const volScalarField& p
) const
{
    return volScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &ThermoType::initializeEnergy,
        p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::calce() const
{
    return volScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &ThermoType::Es,
        this->rho_,
        this->e_,
        this->T_
    );
}

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
        "he",
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
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::Es,
        cells,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        basicBlastThermo::cellSetScalarList(this->e_, cells),
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
Foam::eBlastThermo<BasicThermo, ThermoType>::hs() const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &ThermoType::Hs,
        this->rho_,
        this->e_,
        this->T_
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
        "hs",
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
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::Hs,
        cells,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        basicBlastThermo::cellSetScalarList(this->e_, cells),
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
Foam::eBlastThermo<BasicThermo, ThermoType>::ha() const
{
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &ThermoType::Ha,
        this->rho_,
        this->e_,
        this->T_
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
        "ha",
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
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::Ha,
        cells,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        basicBlastThermo::cellSetScalarList(this->e_, cells),
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
        "hc",
        dimEnergy/dimMass,
        &ThermoType::Hf
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::eBlastThermo<BasicThermo, ThermoType>::flameT() const
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
Foam::eBlastThermo<BasicThermo, ThermoType>::THE() const
{
    return volScalarFieldProperty
    (
        "THE",
        dimTemperature,
        &ThermoType::TRhoE,
        this->T_,
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
        "THE",
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
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::TRhoE,
        cells,
        T,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
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
        "W",
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
    return volScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &ThermoType::kappa,
        this->rho_,
        this->e_,
        this->T_
    );
}


// ************************************************************************* //
