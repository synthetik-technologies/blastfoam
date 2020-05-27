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

#include "detonatingFluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingFluidThermo<Thermo>::detonatingFluidThermo
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
        dict.subDict("reactants"),
        dict.subDict("products"),
        master
    ),
    activation_(activationModel::New(rho.mesh(), dict, name)),
    afterburn_(afterburnModel::New(rho.mesh(), dict, name))
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        Thermo::thermoType1::solid()
     && dict.lookupOrDefault<Switch>("calculateDensity", false)
     && this->rho_.time().value() == this->rho_.time().startTime().value()
    )
    {
        volScalarField rhoInit
        (
            Thermo::blendedVolScalarFieldProperty
            (
                IOobject::groupName("rhoInit", basicThermoModel::name_),
                dimDensity,
                &Thermo::thermoType1::initializeRho,
                &Thermo::thermoType2::initializeRho,
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
Foam::detonatingFluidThermo<Thermo>::~detonatingFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    activation_->solve(stepi, ai, bi);
    afterburn_->solve(stepi, ai, bi);
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    activation_->setODEFields
    (
        nSteps,
        oldIs,
        nOld,
        deltaIs,
        nDelta
    );
    afterburn_->setODEFields
    (
        nSteps,
        oldIs,
        nOld,
        deltaIs,
        nDelta
    );
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::clearODEFields()
{
    activation_->clearODEFields();
    afterburn_->clearODEFields();
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::correct()
{
    if (this->master_)
    {
        this->T_ = this->calcT();
        this->p_ = calcP();
        this->p_.max(small);
    }

    if (this->viscous_)
    {
        this->fluidThermoModel::mu_ = Thermo::blendedVolScalarFieldProperty
            (
                "mu",
                dimDynamicViscosity,
                &Thermo::thermoType1::mu,
                &Thermo::thermoType2::mu,
                this->rho_,
                this->e_,
                this->T_
            );

        this->fluidThermoModel::alpha_ = Thermo::blendedVolScalarFieldProperty
            (
                "alpha",
                dimensionSet(1, -1, -1, 0, 0),
                &Thermo::thermoType1::alphah,
                &Thermo::thermoType2::alphah,
                this->rho_,
                this->e_,
                this->T_
            );
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<Thermo>::speedOfSound() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "speedOfSound",
        dimVelocity,
        &Thermo::thermoType1::speedOfSound,
        &Thermo::thermoType2::speedOfSound,
        this->p_,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<Thermo>::speedOfSound(const label patchi) const
{
    return Thermo::blendedPatchFieldProperty
    (
        &Thermo::thermoType1::speedOfSound,
        &Thermo::thermoType2::speedOfSound,
        patchi,
        this->p_.boundaryField()[patchi],
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<Thermo>::calcP() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "p",
        dimPressure,
        &Thermo::thermoType1::p,
        &Thermo::thermoType2::p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<Thermo>::calce() const
{
    tmp<volScalarField> eInit
    (
        Thermo::volScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &Thermo::thermoType1::initializeEnergy,
            this->p_,
            this->rho_,
            this->e_,
            this->T_
        )
    );

    //- Add detonation energy to initially reacted material
    if (this->rho_.time().value() == this->rho_.time().startTime().value())
    {
        eInit.ref() += activation_->e0()*activation_->lambda();
    }

    return eInit;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->p_.mesh().time().timeName(),
                this->p_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (activation_->ESource() + afterburn_->ESource())*this->rho_
        )
    );
}



// ************************************************************************* //
