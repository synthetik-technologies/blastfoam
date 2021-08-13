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

#include "basicFluidBlastThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidBlastThermo<Thermo>::basicFluidBlastThermo
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
        updateRho(this->phaseFluidBlastThermo::p());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidBlastThermo<Thermo>::~basicFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::basicFluidBlastThermo<Thermo>::updateRho(const volScalarField& p)
{
    volScalarField rhoNew
    (
        Thermo::volScalarFieldProperty
        (
            "rho",
            dimDensity,
            &Thermo::initializeRho,
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
void Foam::basicFluidBlastThermo<Thermo>::correct()
{}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::ESource() const
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
Foam::basicFluidBlastThermo<Thermo>::speedOfSound() const
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
                    Thermo::volScalarFieldProperty
                    (
                        "cSqr",
                        sqr(dimVelocity),
                        &Thermo::thermoType::cSqr,
                        this->phaseFluidBlastThermo::p(),
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
Foam::basicFluidBlastThermo<Thermo>::speedOfSound(const label patchi) const
{
    return sqrt
    (
        max
        (
            Thermo::patchFieldProperty
            (
                &Thermo::thermoType::cSqr,
                patchi,
                this->phaseFluidBlastThermo::p().boundaryField()[patchi],
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
Foam::basicFluidBlastThermo<Thermo>::Gamma() const
{
    return Thermo::volScalarFieldProperty
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
Foam::basicFluidBlastThermo<Thermo>::Gamma(const label patchi) const
{
    return Thermo::patchFieldProperty
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
Foam::basicFluidBlastThermo<Thermo>::Gammai(const label celli) const
{
    return Thermo::thermoType::Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::pRhoT() const
{
    return
        Thermo::volScalarFieldProperty
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
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::pRhoTi(const label celli) const
{
    return Thermo::thermoType::p
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::dpdRhoi(const label celli) const
{
    return Thermo::thermoType::dpdRho
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::dpdei(const label celli) const
{
    return Thermo::thermoType::dpde
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::calce(const volScalarField& p) const
{
    return Thermo::volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &Thermo::initializeEnergy,
        p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::basicFluidBlastThermo<Thermo>::mu() const
{
    return Thermo::volScalarFieldProperty
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
