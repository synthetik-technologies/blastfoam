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

#include "mixtureBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
const Foam::PtrList<ThermoType>&
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::constructSpeciesData
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::volScalarFieldProperty
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
            IOobject::groupName(psiName, this->name()),
            this->rho_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(psi, celli)
    {
        psi[celli] = (this->mixture_[celli].*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(pPsi, facei)
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::volScalarFieldSpecieProperty
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
            IOobject::groupName(psiName, this->name()),
            this->rho_.mesh(),
            psiDim
        )
    );

    const ThermoType& thermo = speciesData_[this->species_[s]];

    volScalarField& psi = tPsi.ref();

    forAll(psi, celli)
    {
        psi[celli] = (thermo.*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(pPsi, facei)
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::patchFieldProperty
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellSetProperty
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::mixtureBlastThermo
(
    const word& name,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
:
    multicomponentBlastThermo<BasicThermo>
    (
        name,
        rho,
        e,
        T,
        dict,
        masterName
    ),
    speciesData_(this->species_.size()),
    mixture_(rho.mesh(), this->Ys_, constructSpeciesData(dict), name)
{}


template<class BasicThermo, class ThermoType>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::mixtureBlastThermo
(
    const HashPtrTable<ThermoType, word, string::hash>& thermoData,
    const word& name,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
:
    multicomponentBlastThermo<BasicThermo>
    (
        thermoData.toc(),
        name,
        rho,
        e,
        T,
        dict,
        masterName
    ),
    speciesData_(this->species_.size()),
    mixture_(rho.mesh(), this->Ys_, constructSpeciesData(dict), name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::~mixtureBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
void Foam::mixtureBlastThermo<BasicThermo, ThermoType>::correct()
{
    this->mixture_.update();
}


template<class BasicThermo, class ThermoType>
void Foam::mixtureBlastThermo<BasicThermo, ThermoType>::updateRho()
{
    volScalarField rhoNew
    (
        this->volScalarFieldProperty
        (
            "rhoNew",
            dimDensity,
            &ThermoType::initializeRho,
            this->p_,
            this->rho_,
            this->e_,
            this->T_
        )
    );
    this->rho_.ref() = rhoNew();

    forAll(this->rho_.boundaryField(), patchi)
    {
        forAll(this->rho_.boundaryField()[patchi], facei)
        {
            this->rho_.boundaryFieldRef()[patchi][facei] =
                rhoNew.boundaryField()[patchi][facei];
        }
    }
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::he
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::Es,
        cells,
        this->cellSetScalarList(this->rho_, cells),
        this->cellSetScalarList(this->e_, cells),
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::he
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hs() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("hs", this->name_),
        dimEnergy/dimMass,
        &ThermoType::Hs,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("hs", this->name_),
        dimEnergy/dimMass,
        &ThermoType::Hs,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::Hs,
        cells,
        this->cellSetScalarList(this->rho_, cells),
        this->cellSetScalarList(this->e_, cells),
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hs
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::ha() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("ha", this->name_),
        dimEnergy/dimMass,
        &ThermoType::Ha,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::ha
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &ThermoType::Ha,
        cells,
        this->cellSetScalarList(this->rho_, cells),
        this->cellSetScalarList(this->e_, cells),
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::ha
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hc() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("hc", basicBlastThermo::name_),
        dimEnergy/dimMass,
        &ThermoType::Hf
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::flameT() const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::THE() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("THE", basicBlastThermo::name_),
        dimTemperature,
        &ThermoType::TRhoE,
        this->T_,
        this->rho_,
        this->e_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("THE", basicBlastThermo::masterName_),
        dimTemperature,
        &ThermoType::TRhoE,
        T,
        this->rho_,
        he
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::THE
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::THE
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::THEi
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return this->mixture_[celli].ThermoType::TRhoE(T, this->rho_[celli], he);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cp() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cp", this->name_),
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cp,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cp
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cpi
(
    const label celli
) const
{
    return this->mixture_[celli].Cp
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cv", this->name_),
        dimEnergy/dimMass/dimTemperature,
        &ThermoType::Cv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cv
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cvi
(
    const label celli
) const
{
    return
        this->mixture_[celli].Cv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cpv() const
{
    return volScalarField::New("Cpv", Cv());
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv(T, patchi);
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::W() const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::W(const label patchi) const
{
    return patchFieldProperty
    (
        &ThermoType::W,
        patchi
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Wi(const label celli) const
{
    return this->mixture_[celli].W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::alpha() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("alpha", basicBlastThermo::name_),
        dimensionSet(1, -1, -1, 0, 0),
        &ThermoType::alphah,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieW(const word& s) const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieWi
(
    const word& s,
    const label celli
) const
{
    return speciesData_[this->species_[s]].W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieGamma
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieGammai
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieCp(const word& s) const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieCpi
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieCv(const word& s) const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieCvi
(
    const word& s,
    const label celli
) const
{
    return this->speciesData_[this->species_[s]].Cv
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}

template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieCpByCv
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieHf(const word& s) const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::specieFlameT
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
