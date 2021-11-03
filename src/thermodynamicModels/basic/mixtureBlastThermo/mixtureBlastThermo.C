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
#include "blastThermo.H"

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
            IOobject::groupName(psiName, BasicThermo::phaseName()),
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
    const label speciei,
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
            IOobject::groupName(psiName, BasicThermo::phaseName()),
            this->rho_.mesh(),
            psiDim
        )
    );

    const ThermoType& thermo = speciesData_[speciei];

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
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    multicomponentBlastThermo
    (
        mesh,
        dict,
        phaseName,
        masterName
    ),
    BasicThermo
    (
        mesh,
        dict,
        phaseName,
        masterName
    ),
    speciesData_(this->species_.size()),
    mixture_(mesh, this->Y_, constructSpeciesData(dict), phaseName)
{}


template<class BasicThermo, class ThermoType>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::mixtureBlastThermo
(
    const HashPtrTable<ThermoType, word, string::hash>& thermoData,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    multicomponentBlastThermo
    (
        thermoData.toc(),
        mesh,
        dict,
        phaseName,
        masterName
    ),
    BasicThermo
    (
        thermoData.toc(),
        mesh,
        dict,
        phaseName,
        masterName
    ),
    speciesData_(this->species_.size()),
    mixture_(mesh, this->Y_, constructSpeciesData(dict), phaseName)
{}


template<class BasicThermo, class ThermoType>
void Foam::mixtureBlastThermo<BasicThermo, ThermoType>::initializeModels()
{
    BasicThermo::initializeModels();
    multicomponentBlastThermo::initializeModels();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::~mixtureBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasicThermo, class ThermoType>
void Foam::mixtureBlastThermo<BasicThermo, ThermoType>::updateMixture()
{
    this->mixture_.updateMixture();
}


template<class BasicThermo, class ThermoType>
void Foam::mixtureBlastThermo<BasicThermo, ThermoType>::solve()
{
    multicomponentBlastThermo::solve();
}



template<class BasicThermo, class ThermoType>
void Foam::mixtureBlastThermo<BasicThermo, ThermoType>::postUpdate()
{
    multicomponentBlastThermo::postUpdate();
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
        blastThermo::cellSetScalarList(this->rho_, cells),
        blastThermo::cellSetScalarList(this->e_, cells),
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
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellHE
(
    const scalar T,
    const label celli
) const
{
    return this->mixture_[celli].ThermoType::Es
    (
        this->rho_[celli],
        this->e_[celli],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::patchFaceHE
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    return this->mixture_.boundary(patchi, facei).ThermoType::Es
    (
        this->rho_.boundaryField()[patchi][facei],
        this->e_.boundaryField()[patchi][facei],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hs() const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::hs
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
        blastThermo::cellSetScalarList(this->rho_, cells),
        blastThermo::cellSetScalarList(this->e_, cells),
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::ha
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
        blastThermo::cellSetScalarList(this->rho_, cells),
        blastThermo::cellSetScalarList(this->e_, cells),
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
        "hc",
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
        "flameT",
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::THE
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellTHE
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return this->mixture_[celli].ThermoType::TRhoE(T, this->rho_[celli], he);
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellCp
(
    const scalar T,
    const label celli
) const
{
    return this->mixture_[celli].Cp
    (
        this->rho_[celli],
        this->e_[celli],
        T
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellCv
(
    const scalar T,
    const label celli
) const
{
    return
        this->mixture_[celli].Cv
        (
            this->rho_[celli],
            this->e_[celli],
            T
        );
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
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellCpv
(
    const scalar T,
    const label celli
) const
{
    return this->mixture_[celli].ThermoType::Cv
    (
        this->rho_[celli],
        this->e_[celli],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::patchFaceCpv
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    return this->mixture_.boundary(patchi, facei).ThermoType::Cv
    (
        this->rho_.boundaryField()[patchi][facei],
        this->e_.boundaryField()[patchi][facei],
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::W() const
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::cellW(const label celli) const
{
    return this->mixture_[celli].W();
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::calcKappa() const
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


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Wi
(
    const label speciei
) const
{
    return speciesData_[speciei].W();
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Hf
(
    const label speciei
) const
{
    return speciesData_[speciei].Hf();
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::flameT
(
    const label speciei
) const
{
    return speciesData_[speciei].flameT();
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::rho
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return speciesData_[speciei].rhoPT
    (
        p,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::rho
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldSpecieProperty
    (
        speciei,
        "rho",
        dimDensity,
        &ThermoType::rhoPT,
        p,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cp
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Cp
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldSpecieProperty
    (
        speciei,
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
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::HE
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::HE
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::HE
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldSpecieProperty
    (
        speciei,
        "HE",
        dimEnergy/dimMass,
        &ThermoType::Es,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Hs
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Hs
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Hs
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldSpecieProperty
    (
        speciei,
        "Hs",
        dimEnergy/dimMass,
        &ThermoType::Hs,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Ha
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Ha
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::Ha
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldSpecieProperty
    (
        speciei,
        "Ha",
        dimEnergy/dimMass,
        &ThermoType::Ha,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::mu
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::kappa
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::mixtureBlastThermo<BasicThermo, ThermoType>::kappa
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldSpecieProperty
    (
        speciei,
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &ThermoType::kappa,
        this->rho_,
        this->e_,
        T
    );
}

// ************************************************************************* //
