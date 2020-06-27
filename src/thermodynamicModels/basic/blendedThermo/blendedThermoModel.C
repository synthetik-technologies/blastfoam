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

#include "blendedThermoModel.H"

// * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::volScalarFieldProperty
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


template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::cellSetProperty
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


template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::patchFieldProperty
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


template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method1, class Method2, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::blendedVolScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    Method1 psiMethod1,
    Method2 psiMethod2,
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
        scalar x = this->xi(celli);
        if (x < small)
        {
            psi[celli] = (this->*psiMethod1)(args[celli] ...);
        }
        else if ((1.0 - x) < small)
        {
            psi[celli] = (this->*psiMethod2)(args[celli] ...);
        }
        else
        {
            psi[celli] =
                (this->*psiMethod2)(args[celli] ...)*x
              + (this->*psiMethod1)(args[celli] ...)*(1.0 - x);
        }
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];
        tmp<scalarField> xp(this->x(patchi));

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            const scalar& x = xp()[facei];
            if (x < small)
            {
                pPsi[facei] =
                    (this->*psiMethod1)(args.boundaryField()[patchi][facei] ...);
            }
            else if ((1.0 - x) < small)
            {
                pPsi[facei] =
                    (this->*psiMethod2)(args.boundaryField()[patchi][facei] ...);
            }
            else
            {
                pPsi[facei] =
                    (this->*psiMethod2)(args.boundaryField()[patchi][facei] ...)*x
                  + (this->*psiMethod1)
                    (
                        args.boundaryField()[patchi][facei] ...
                    )*(1.0 - x);
            }
        }
    }

    return tPsi;
}


template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method1, class Method2, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::blendedCellSetProperty
(
    Method1 psiMethod1,
    Method2 psiMethod2,
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
        scalar x = this->xi(cells[celli]);
        if (x < small)
        {
            psi[celli] = (this->*psiMethod1)(args[celli] ...);
        }
        else if ((1.0 - x) < small)
        {
            psi[celli] = (this->*psiMethod2)(args[celli] ...);
        }
        else
        {
            psi[celli] =
                (this->*psiMethod2)(args[celli] ...)*x
              + (this->*psiMethod1)(args[celli] ...)*(1.0 - x);
        }
    }

    return tPsi;
}


template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method1, class Method2, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::blendedPatchFieldProperty
(
    Method1 psiMethod1,
    Method2 psiMethod2,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();
    tmp<scalarField> xtmp(this->x(patchi));
    const scalarField& x(xtmp());

    forAll(this->p_.boundaryField()[patchi], facei)
    {
        if (x[facei] < small)
        {
             psi[facei] = (this->*psiMethod1)(args[facei] ...);
        }
        else if ((1.0 - x[facei]) < small)
        {
            psi[facei] = (this->*psiMethod2)(args[facei] ...);
        }
        else
        {
            psi[facei] =
                (this->*psiMethod2)(args[facei] ...)*x[facei]
              + (this->*psiMethod1)(args[facei] ...)*(1.0 - x[facei]);
        }
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class Thermo1, class Thermo2>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::blendedThermoModel
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const dictionary& dict1,
    const dictionary& dict2,
    const bool master
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
        master
    ),
    Thermo1(dict1),
    Thermo2(dict2)
{}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::blendedThermoModel
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& dict1,
    const dictionary& dict2,
    const bool master
)
:
    BasicThermo
    (
        name,
        mesh,
        dict,
        master
    ),
    Thermo1(dict1),
    Thermo2(dict2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class Thermo1, class Thermo2>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::~blendedThermoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::calcT() const
{
    return blendedVolScalarFieldProperty
    (
        IOobject::groupName("T", basicThermoModel::name_),
        dimTemperature,
        &Thermo1::TRhoE,
        &Thermo2::TRhoE,
        this->T_,
        this->rho_,
        this->e_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::TRhoE
(
    const scalarField& T,
    const scalarField& e,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::TRhoE,
        &Thermo2::TRhoE,
        patchi,
        T,
        this->rho_.boundaryField()[patchi],
        e
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::TRhoEi
(
    const scalar& T,
    const scalar& e,
    const label celli
) const
{
    return
        Thermo1::TRhoE
        (
            this->T_[celli],
            this->rho_[celli],
            this->e_[celli]
        )*this->xi(celli)
      + Thermo2::TRhoE
        (
            this->T_[celli],
            this->rho_[celli],
            this->e_[celli]
        )*(1.0 - this->xi(celli));
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::E() const
{
    return blendedVolScalarFieldProperty
    (
        IOobject::groupName("e", basicThermoModel::name_),
        dimEnergy/dimMass,
        &Thermo1::Es,
        &Thermo2::Es,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::e
(
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    return blendedVolScalarFieldProperty
    (
        IOobject::groupName("e", basicThermoModel::name_),
        dimEnergy/dimMass,
        &Thermo1::Es,
        &Thermo2::Es,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Es,
        &Thermo2::Es,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    return blendedCellSetProperty
    (
        &Thermo1::Es,
        &Thermo2::Es,
        faceCells,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::W() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            "W",
            this->p_.mesh(),
            dimMass/dimMoles
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(this->p_, celli)
    {
        scalar x = this->xi(celli);
        if (x < small)
        {
            psi[celli] = Thermo1::W();
        }
        else if ((1.0 - x) < small)
        {
            psi[celli] = Thermo2::W();
        }
        else
        {
            psi[celli] = Thermo1::W()*x + Thermo2::W()*(1.0 - x);
        }
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];
        tmp<scalarField> xp(this->x(patchi));

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            const scalar& x = xp()[facei];
            if (x < small)
            {
                pPsi[facei] = Thermo1::W();
            }
            else if ((1.0 - x) < small)
            {
                pPsi[facei] = Thermo2::W();
            }
            else
            {
                pPsi[facei] = Thermo1::W()*x + Thermo2::W()*(1.0 - x);
            }
        }
    }

    return tPsi;
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::W(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();
    const scalarField x(this->x(patchi));

    forAll(this->p_.boundaryField()[patchi], facei)
    {
        if (x[facei] < small)
        {
             psi[facei] = Thermo1::W();
        }
        else if ((1.0 - x[facei]) < small)
        {
            psi[facei] = Thermo2::W();
        }
        else
        {
            psi[facei] =Thermo1::W()*x[facei] + Thermo2::W()*(1.0 - x[facei]);
        }
    }

    return tPsi;
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Wi(const label celli) const
{
    return this->xi(celli)*Thermo1::W() + (1.0 - this->xi(celli))*Thermo2::W();
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Gamma() const
{
    return blendedVolScalarFieldProperty
    (
        "Gamma",
        dimless,
        &Thermo1::Gamma,
        &Thermo2::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Gamma(const label patchi) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Gamma,
        &Thermo2::Gamma,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cp() const
{
    return blendedVolScalarFieldProperty
    (
        "Cp",
        dimEnergy/dimMass/dimTemperature,
        &Thermo1::Cp,
        &Thermo2::Cp,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cp(const label patchi) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Cp,
        &Thermo2::Cp,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Cp,
        &Thermo2::Cp,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cpi(const label celli) const
{
    return
        Thermo1::Cp
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        )*this->xi(celli)
      + Thermo2::Cp
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        )*(1.0 - this->xi(celli));
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cv() const
{
    return blendedVolScalarFieldProperty
    (
        "Cv",
        dimEnergy/dimMass/dimTemperature,
        &Thermo1::Cv,
        &Thermo2::Cv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cv(const label patchi) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Cv,
        &Thermo2::Cv,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Cv,
        &Thermo2::Cv,
        patchi,
        rho,
        e,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::Cvi(const label celli) const
{
    return
        Thermo1::Cv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        )*this->xi(celli)
      + Thermo2::Cv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        )*(1.0 - this->xi(celli));
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::CpByCv() const
{
    return blendedVolScalarFieldProperty
    (
        "CpByCv",
        dimless,
        &Thermo1::CpByCv,
        &Thermo2::CpByCv,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::CpByCv(const label patchi) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::CpByCv,
        &Thermo2::CpByCv,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedThermoModel<BasicThermo, Thermo1, Thermo2>::CpByCv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::CpByCv,
        &Thermo2::CpByCv,
        patchi,
        rho,
        e,
        T
    );
}

// ************************************************************************* //
