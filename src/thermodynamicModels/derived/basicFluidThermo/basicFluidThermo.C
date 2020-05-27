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

#include "basicFluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidThermo<Thermo>::basicFluidThermo
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
    Thermo
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master
    )
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        dict.lookupOrDefault<Switch>("calculateDensity", false)
     && this->rho_.time().value() == this->rho_.time().startTime().value()
    )
    {
        volScalarField rhoInit
        (
            Thermo::volScalarFieldProperty
            (
                "rhoInit",
                dimDensity,
                &Thermo::initializeRho,
                this->p_,
                this->rho_,
                this->e_,
                this->T_
            )
        );
        this->rho_ = rhoInit;
        forAll(this->rho_.boundaryField(), patchi)
        {
            forAll(this->rho_.boundaryField()[patchi], facei)
            {
                this->rho_.boundaryFieldRef()[patchi][facei] =
                    rhoInit.boundaryField()[patchi][facei];
            }
        }
    }

    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidThermo<Thermo>::~basicFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::basicFluidThermo<Thermo>::correct()
{
    if (this->master_)
    {
        this->T_ = this->calcT();
        this->p_ = calcP();
        this->p_.max(small);
    }

    if (this->viscous_)
    {
        this->fluidThermoModel::mu_ = Thermo::volScalarFieldProperty
        (
            "thermo:mu",
            dimDynamicViscosity,
            &Thermo::thermoType::mu,
            this->rho_,
            this->e_,
            this->T_
        );

        this->fluidThermoModel::alpha_ = Thermo::volScalarFieldProperty
        (
            "alphah",
            dimensionSet(1, -1, -1, 0, 0),
            &Thermo::thermoType::alphah,
            this->rho_,
            this->e_,
            this->T_
        );
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->p_.mesh().time().timeName(),
                this->p_.mesh()
            ),
            this->p_.mesh(),
            dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::speedOfSound() const
{
    return Thermo::volScalarFieldProperty
    (
        "speedOfSound",
        dimVelocity,
        &Thermo::thermoType::speedOfSound,
        this->p_,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::basicFluidThermo<Thermo>::speedOfSound(const label patchi) const
{
    return Thermo::patchFieldProperty
    (
        &Thermo::thermoType::speedOfSound,
        patchi,
        this->p_.boundaryField()[patchi],
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::calcP() const
{
    return Thermo::volScalarFieldProperty
    (
        "p",
        dimPressure,
        &Thermo::thermoType::p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidThermo<Thermo>::calce() const
{
    return Thermo::volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &Thermo::initializeEnergy,
        this->p_,
        this->rho_,
        this->e_,
        this->T_
    );
}

// ************************************************************************* //
