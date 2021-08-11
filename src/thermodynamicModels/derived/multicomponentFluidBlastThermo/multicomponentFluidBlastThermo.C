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
    const word& name,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
:
    Thermo
    (
        name,
        rho,
        e,
        T,
        dict,
        masterName
    )
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() == 0
     || (
            dict.lookupOrDefault<Switch>("calculateDensity", false)
         && this->rho_.time().timeIndex() == 0
        )
    )
    {
        updateRho(phaseFluidBlastThermo::p());
    }
}


template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::multicomponentFluidBlastThermo
(
    const HashPtrTable<Thermo, word, string::hash>& thermoData,
    const word& name,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
:
    Thermo
    (
        thermoData,
        name,
        rho,
        e,
        T,
        dict,
        masterName
    )
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() == 0
     || (
            dict.lookupOrDefault<Switch>("calculateDensity", false)
         && this->rho_.time().timeIndex() == 0
        )
    )
    {
        updateRho();
    }

    this->mu_ =
        this->volScalarFieldProperty
        (
            IOobject::groupName("mu", name),
            dimDynamicViscosity,
            &Thermo::thermoType::mu,
            this->rho_,
            this->e_,
            this->T_
        );

    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::~multicomponentFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentFluidBlastThermo<Thermo>::updateRho(const volScalarField& p)
{
    volScalarField rhoNew
    (
        Thermo::volScalarFieldProperty
        (
            IOobject::groupName("rhoNew", this->name()),
            dimDensity,
            &Thermo::thermoType::initializeRho,
            p,
            this->rho_,
            this->e_,
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
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::speedOfSound() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("speedOfSound", this->name()),
            sqrt
            (
                max
                (
                    this->volScalarFieldProperty
                    (
                        "cSqr",
                        sqr(dimVelocity),
                        &Thermo::thermoType::cSqr,
                        phaseFluidBlastThermo::p(),
                        this->rho_,
                        this->e_,
                        this->T_
                    ),
                    dimensionedScalar(sqr(dimVelocity), small)
                )
            )
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::speedOfSound(const label patchi) const
{
    return sqrt
    (
        max
        (
            this->patchFieldProperty
            (
                &Thermo::thermoType::cSqr,
                patchi,
                phaseFluidBlastThermo::p().boundaryField()[patchi],
                this->rho_.boundaryField()[patchi],
                this->e_.boundaryField()[patchi],
                this->T_.boundaryField()[patchi]
            ),
            small
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::Gamma() const
{
    return this->volScalarFieldProperty
    (
        "Gamma",
        dimless,
        &Thermo::thermoType::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::Gamma(const label patchi) const
{
    return this->patchFieldProperty
    (
        &Thermo::thermoType::Gamma,
        patchi,
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::Gammai(const label celli) const
{
    return this->mixture_[celli].Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::pRhoT() const
{
    return
        this->volScalarFieldProperty
        (
            IOobject::groupName("p", this->name_),
            dimPressure,
            &Thermo::thermoType::p,
            this->rho_,
            this->e_,
            this->T_
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
Foam::multicomponentFluidBlastThermo<Thermo>::mu() const
{
    return this->volScalarFieldProperty
    (
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        &Thermo::thermoType::mu,
        this->rho_,
        this->e_,
        this->T_
    );
}

// ************************************************************************* //
