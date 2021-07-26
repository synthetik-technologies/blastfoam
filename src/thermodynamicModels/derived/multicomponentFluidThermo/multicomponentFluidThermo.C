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

#include "multicomponentFluidThermo.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multicomponentFluidThermo<ThermoType>::multicomponentFluidThermo
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    mixtureThermoModel<multicomponentFluidThermoModel, ThermoType>
    (
        name,
        mesh,
        dict,
        master,
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
            &ThermoType::mu,
            this->rho_,
            this->e_,
            this->T_
        );

    this->initialize();
}


template<class ThermoType>
Foam::multicomponentFluidThermo<ThermoType>::multicomponentFluidThermo
(
    const HashPtrTable<ThermoType, word, string::hash>& thermoData,
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    mixtureThermoModel<multicomponentFluidThermoModel, ThermoType>
    (
        thermoData,
        name,
        mesh,
        dict,
        master,
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
            &ThermoType::mu,
            this->rho_,
            this->e_,
            this->T_
        );

    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multicomponentFluidThermo<ThermoType>::~multicomponentFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::multicomponentFluidThermo<ThermoType>::updateRho()
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

template<class ThermoType>
void Foam::multicomponentFluidThermo<ThermoType>::correct()
{
    this->mixture_.update();

    this->T_ = this->calcT();
    this->T_.correctBoundaryConditions();

    this->p_ = fluidThermoModel::calcP();
    this->p_.max(small);
    this->p_.correctBoundaryConditions();

    if (this->viscous_)
    {
        this->mu_ =
            this->volScalarFieldProperty
            (
                "thermo:mu",
                dimDynamicViscosity,
                &ThermoType::mu,
                this->rho_,
                this->e_,
                this->T_
            );

        this->alpha_ =
            this->volScalarFieldProperty
            (
                "alphah",
                dimensionSet(1, -1, -1, 0, 0),
                &ThermoType::alphah,
                this->rho_,
                this->e_,
                this->T_
            );
    }
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<ThermoType>::ESource() const
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


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<ThermoType>::speedOfSound() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("speedOfSound", this->group()),
            sqrt
            (
                max
                (
                    this->volScalarFieldProperty
                    (
                        "cSqr",
                        sqr(dimVelocity),
                        &ThermoType::thermoType::cSqr,
                        this->p_,
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


template<class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<ThermoType>::speedOfSound(const label patchi) const
{
    return sqrt
    (
        max
        (
            this->patchFieldProperty
            (
                &ThermoType::thermoType::cSqr,
                patchi,
                this->p_.boundaryField()[patchi],
                this->rho_.boundaryField()[patchi],
                this->e_.boundaryField()[patchi],
                this->T_.boundaryField()[patchi]
            ),
            small
        )
    );
}


template<class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<ThermoType>::calcP(const label patchi) const
{
    return
        this->patchFieldProperty
        (
            &ThermoType::p,
            patchi,
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi]
        );
}


template<class ThermoType>
Foam::scalar Foam::multicomponentFluidThermo<ThermoType>::calcPi(const label celli) const
{
    return
        this->mixture_[celli].p
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<ThermoType>::calce() const
{
    return this->volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &ThermoType::initializeEnergy,
        this->p_,
        this->rho_,
        this->e_,
        this->T_
    );
}

// ************************************************************************* //
