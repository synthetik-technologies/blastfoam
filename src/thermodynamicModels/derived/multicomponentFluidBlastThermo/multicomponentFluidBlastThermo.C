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

#include "multicomponentFluidBlastThermo.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::multicomponentFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    Thermo
    (
        mesh,
        dict,
        phaseName,
        masterName
    ),
    localMixture_(this->mixture_)
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() <= 0
     || (
            dict.lookupOrDefault<Switch>("calculateDensity", false)
         && !this->rho_.time().restart()
        )
    )
    {
        updateRho(Thermo::baseThermo::p());
    }
    this->initializeFields();
}


template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::multicomponentFluidBlastThermo
(
    const HashPtrTable<Thermo, word, string::hash>& thermoData,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    Thermo
    (
        thermoData,
        mesh,
        dict,
        phaseName,
        masterName
    ),
    localMixture_(this->mixture_)
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() <= 0
     || (
            dict.lookupOrDefault<Switch>("calculateDensity", false)
         && !this->rho_.time().restart()
        )
    )
    {
        updateRho(Thermo::baseThermo::p());
    }
    this->initializeFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::~multicomponentFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentFluidBlastThermo<Thermo>::updateRho
(
    const volScalarField& p
)
{
    volScalarField rhoNew
    (
        Thermo::volScalarFieldProperty
        (
            "rhoNew",
            dimDensity,
            &Thermo::thermoType::rhoPT,
            p,
            this->T_
        )
    );
    this->rho_ = rhoNew;
    forAll(this->rho_.boundaryField(), patchi)
    {
        forAll(this->rho_.boundaryField()[patchi], facei)
        {
            this->rho_.boundaryFieldRef()[patchi][facei] =
                rhoNew.boundaryField()[patchi][facei];
        }
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh()
            ),
            this->rho_.mesh(),
            dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::cellpRhoT(const label celli) const
{
    return this->mixture_[celli].p
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    return this->mixture_[celli].Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}

template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    return this->mixture_[celli].dpdRho
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::celldpde(const label celli) const
{
    return this->mixture_[celli].dpde
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::calce
(
    const volScalarField& p
) const
{
    return this->volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &Thermo::thermoType::initializeEnergy,
        p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::mu
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return this->volScalarFieldSpecieProperty
    (
        speciei,
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        &Thermo::thermoType::mu,
        this->rho_,
        this->e_,
        T
    );
}

// ************************************************************************* //
