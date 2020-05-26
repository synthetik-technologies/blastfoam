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

#include "SolidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Thermo>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField> Foam::SolidThermoModel<Thermo>::volScalarFieldProperty
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
Foam::tmp<Foam::scalarField> Foam::SolidThermoModel<Thermo>::cellSetProperty
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
Foam::tmp<Foam::scalarField> Foam::SolidThermoModel<Thermo>::patchFieldProperty
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
Foam::SolidThermoModel<Thermo>::SolidThermoModel
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
    solidThermoModel
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
    volScalarField& rhoRef(solidThermoModel::rho_);
    forAll(rhoRef, celli)
    {
        rhoRef[celli] = Thermo::rho();
    }
    forAll(rhoRef.boundaryField(), patchi)
    {
        forAll(rhoRef.boundaryField()[patchi], facei)
        {
            rhoRef.boundaryFieldRef()[patchi][facei] = Thermo::rho();
        }
    }
    e_ = E();
    T_ = calcT();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::SolidThermoModel<Thermo>::~SolidThermoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::SolidThermoModel<Thermo>::correct()
{
    if (master_)
    {
        T_ = calcT();
    }
}

template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::calcT() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("T", basicThermoModel::name_),
        dimTemperature,
        &Thermo::TRhoE,
        T_,
        rho_,
        e_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::E() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("e", basicThermoModel::name_),
        dimEnergy/dimMass,
        &Thermo::Es,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoModel<Thermo>::E
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &Thermo::Es,
        patchi,
        rho,
        e,
        T
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoModel<Thermo>::E
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &Thermo::Es,
        faceCells,
        rho,
        e,
        T
    );
}

template<class Thermo>
Foam::tmp<Foam::volVectorField>
Foam::SolidThermoModel<Thermo>::Kappa() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volVectorField> tKappa
    (
        volVectorField::New
        (
            "Kappa",
            mesh,
            dimEnergy/dimTime/dimLength/dimTemperature
        )
    );

    volVectorField& Kappa = tKappa.ref();
    vectorField& KappaCells = Kappa.primitiveFieldRef();
    const scalarField& rhoCells = solidThermoModel::rho_;
    const scalarField& eCells = this->e_;
    const scalarField& TCells = this->T_;

    forAll(KappaCells, celli)
    {
        Kappa[celli] = Thermo::Kappa(rhoCells[celli], eCells[celli], TCells[celli]);
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        vectorField& Kappap = KappaBf[patchi];
        const scalarField& pRho = solidThermoModel::rho_.boundaryField()[patchi];
        const scalarField& pe = this->e_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] = Thermo::Kappa(pRho[facei], pe[facei], pT[facei]);
        }
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::vectorField>
Foam::SolidThermoModel<Thermo>::Kappa(const label patchi) const
{
    const scalarField& pRho = solidThermoModel::rho_.boundaryField()[patchi];
    const scalarField& pe = this->e_.boundaryField()[patchi];
    const scalarField& pT = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pe.size()));

    vectorField& Kappap = tKappa.ref();

    forAll(pe, facei)
    {
        Kappap[facei] = Thermo::Kappa(pRho[facei], pe[facei], pT[facei]);
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::Gamma() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Gamma", basicThermoModel::name_),
        dimless,
        &Thermo::Gamma,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoModel<Thermo>::Gamma(const label patchi) const
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
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::W() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("W", basicThermoModel::name_),
        dimMass/dimMoles,
        &Thermo::W
    );
}


template<class Thermo>
Foam::scalar
Foam::SolidThermoModel<Thermo>::Wi(const label celli) const
{
    return Thermo::W();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::Cp() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cp", basicThermoModel::name_),
        dimEnergy/dimMass/dimTemperature,
        &Thermo::Cp,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoModel<Thermo>::Cp(const label patchi) const
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
Foam::scalar
Foam::SolidThermoModel<Thermo>::Cpi(const label celli) const
{
    return Thermo::Cp(rho_[celli], e_[celli], T_[celli]);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::Cv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cv", basicThermoModel::name_),
        dimEnergy/dimMass/dimTemperature,
        &Thermo::Cv,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoModel<Thermo>::Cv(const label patchi) const
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
Foam::SolidThermoModel<Thermo>::Cv
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
Foam::scalar
Foam::SolidThermoModel<Thermo>::Cvi(const label celli) const
{
    return Thermo::Cv(rho_[celli], e_[celli], T_[celli]);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::SolidThermoModel<Thermo>::CpByCv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("CpByCv", basicThermoModel::name_),
        dimless,
        &Thermo::CpByCv,
        rho_,
        e_,
        T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::SolidThermoModel<Thermo>::CpByCv(const label patchi) const
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
