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

#include "blendedBlastThermo.H"

// * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::volScalarFieldProperty
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
            IOobject::groupName(psiName, this->phaseName()),
            this->rho_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(psi, celli)
    {
        psi[celli] = (this->*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(pPsi, facei)
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::cellSetProperty
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::patchFieldProperty
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
            (this->*psiMethod)(args[facei] ...);
    }

    return tPsi;
}


template<class BasicThermo, class Thermo1, class Thermo2>
template<class Method1, class Method2, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::blendedVolScalarFieldProperty
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
            IOobject::groupName(psiName, this->phaseName()),
            this->rho_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(psi, celli)
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

        forAll(pPsi, facei)
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::blendedCellSetProperty
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::blendedPatchFieldProperty
(
    Method1 psiMethod1,
    Method2 psiMethod2,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->rho_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();
    tmp<scalarField> xtmp(this->x(patchi));
    const scalarField& x(xtmp());

    forAll(psi, facei)
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::blendedBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& dict1,
    const dictionary& dict2,
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
    Thermo1(dict1),
    Thermo2(dict2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class Thermo1, class Thermo2>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::~blendedBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class Thermo1, class Thermo2>
Foam::word
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::thermoName() const
{
    return Thermo1::typeName() + "/" + Thermo2::typeName();
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return blendedVolScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &Thermo1::Es,
        &Thermo2::Es,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return blendedCellSetProperty
    (
        &Thermo1::Es,
        &Thermo2::Es,
        cells,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        basicBlastThermo::cellSetScalarList(this->e_, cells),
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::he
(
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Es,
        &Thermo2::Es,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::hs() const
{
    return blendedVolScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &Thermo1::Hs,
        &Thermo2::Hs,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return blendedVolScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &Thermo1::Hs,
        &Thermo2::Hs,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return blendedCellSetProperty
    (
        &Thermo1::Hs,
        &Thermo2::Hs,
        cells,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        basicBlastThermo::cellSetScalarList(this->e_, cells),
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Hs,
        &Thermo2::Hs,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::ha() const
{
    return blendedVolScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &Thermo1::Ha,
        &Thermo2::Ha,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return blendedVolScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &Thermo1::Ha,
        &Thermo2::Ha,
        this->rho_,
        this->e_,
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return blendedCellSetProperty
    (
        &Thermo1::Ha,
        &Thermo2::Ha,
        cells,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        basicBlastThermo::cellSetScalarList(this->e_, cells),
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Ha,
        &Thermo2::Ha,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::hc() const
{
    return blendedVolScalarFieldProperty
    (
        "hc",
        dimEnergy/dimMass,
        &Thermo1::Hf,
        &Thermo2::Hf
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::flameT() const
{
    return blendedVolScalarFieldProperty
    (
        "flameT",
        dimTemperature,
        &Thermo1::flameT,
        &Thermo2::flameT
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::THE() const
{
    return blendedVolScalarFieldProperty
    (
        "THE",
        dimTemperature,
        &Thermo1::TRhoE,
        &Thermo2::TRhoE,
        this->T_,
        this->rho_,
        this->e_
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return blendedVolScalarFieldProperty
    (
        "THE",
        dimTemperature,
        &Thermo1::TRhoE,
        &Thermo2::TRhoE,
        T,
        this->rho_,
        he
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    return blendedCellSetProperty
    (
        &Thermo1::TRhoE,
        &Thermo2::TRhoE,
        cells,
        T,
        basicBlastThermo::cellSetScalarList(this->rho_, cells),
        he
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::TRhoE,
        &Thermo2::Ha,
        patchi,
        T,
        this->rho_.boundaryField()[patchi],
        he
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::cellTHE
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    const scalar& x = this->xi(celli);
    if (x < small)
    {
        return Thermo1::TRhoE(T, this->rho_[celli], he);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo2::TRhoE(T, this->rho_[celli], he);
    }

    return
        Thermo2::TRhoE(T, this->rho_[celli], he)*x
      + Thermo1::TRhoE(T, this->rho_[celli], he)*(1.0 - x);
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::Cp() const
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Cp,
        &Thermo2::Cp,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::cellCp
(
    const scalar T,
    const label celli
) const
{
    const scalar& x = this->xi(celli);
    if (x < small)
    {
        return Thermo1::Cp(this->rho_[celli], this->e_[celli], T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo2::Cp(this->rho_[celli], this->e_[celli], T);
    }

    return
        Thermo2::Cp(this->rho_[celli], this->e_[celli], T)*x
      + Thermo1::Cp(this->rho_[celli], this->e_[celli], T)*(1.0 - x);
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::Cv() const
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
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::Cv,
        &Thermo2::Cv,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::cellCv
(
    const scalar T,
    const label celli
) const
{
    const scalar& x = this->xi(celli);
    if (x < small)
    {
        return Thermo1::Cv(this->rho_[celli], this->e_[celli], T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo2::Cv(this->rho_[celli], this->e_[celli], T);
    }

    return
        Thermo2::Cv(this->rho_[celli], this->e_[celli], T)*x
      + Thermo1::Cv(this->rho_[celli], this->e_[celli], T)*(1.0 - x);
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::Cpv() const
{
    return volScalarField::New("Cpv", Cv());
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv(T, patchi);
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::W() const
{
    return blendedVolScalarFieldProperty
    (
        "W",
        dimMass/dimMoles,
        &Thermo1::W,
        &Thermo2::W
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::scalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::W(const label patchi) const
{
    return blendedPatchFieldProperty
    (
        &Thermo1::W,
        &Thermo2::W,
        patchi
    );
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::scalar
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::cellW
(
    const label celli
) const
{
    return
        this->xi(celli)*Thermo2::W()
      + (1.0 - this->xi(celli))*Thermo1::W();
}


template<class BasicThermo, class Thermo1, class Thermo2>
Foam::tmp<Foam::volScalarField>
Foam::blendedBlastThermo<BasicThermo, Thermo1, Thermo2>::kappa() const
{
    return blendedVolScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &Thermo1::kappa,
        &Thermo2::kappa,
        this->rho_,
        this->e_,
        this->T_
    );
}


// ************************************************************************* //
