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
    const bool master,
    const word& masterName
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
        master,
        masterName
    )
{
    this->fluidBlastThermo::mu_ =
        max
        (
            Thermo::volScalarFieldProperty
            (
                "mu",
                dimDynamicViscosity,
                &Thermo::thermoType1::mu,
                this->rho_,
                this->e_,
                this->T_
            ),
            Thermo::volScalarFieldProperty
            (
                "mu",
                dimDynamicViscosity,
                &Thermo::thermoType2::mu,
                this->rho_,
                this->e_,
                this->T_
            )
        );
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::initializeModels()
{
    activation_ = activationModel::New
    (
        this->mesh_,
        this->thermoDict_,
        this->name_
    );
    afterburn_ = afterburnModel::New
    (
        this->mesh_,
        this->thermoDict_,
        this->name_
    );
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() == 0
     || (
            Thermo::thermoType1::solid()
         && this->thermoDict_.template lookupOrDefault<Switch>
            (
                "calculateDensity",
                false
            )
         && this->rho_.time().timeIndex() == 0
        )
    )
    {
        updateRho();
    }

    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingFluidThermo<Thermo>::~detonatingFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::update()
{
    activation_->update();
    afterburn_->update();
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::solve()
{
    activation_->solve();
    afterburn_->solve();
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::postUpdate()
{
    activation_->postUpdate();
    afterburn_->postUpdate();
}


template<class Thermo>
void Foam::detonatingFluidThermo<Thermo>::updateRho()
{
    volScalarField rhoNew
    (
        Thermo::blendedVolScalarFieldProperty
        (
            IOobject::groupName("rhoNew", basicBlastThermo::name_),
            dimDensity,
            &Thermo::thermoType1::initializeRho,
            &Thermo::thermoType2::initializeRho,
            this->p_,
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
void Foam::detonatingFluidThermo<Thermo>::correct()
{
    if (this->master_)
    {
        this->T_ = this->calcT();
        this->T_.correctBoundaryConditions();
        this->p_ = fluidBlastThermo::calcP();
        this->p_.correctBoundaryConditions();
    }

    if (this->viscous_)
    {
        this->fluidBlastThermo::mu_ = Thermo::blendedVolScalarFieldProperty
            (
                "mu",
                dimDynamicViscosity,
                &Thermo::thermoType1::mu,
                &Thermo::thermoType2::mu,
                this->rho_,
                this->e_,
                this->T_
            );

        this->fluidBlastThermo::alpha_ = Thermo::blendedVolScalarFieldProperty
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
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("speedOfSound", this->group()),
            sqrt
            (
                max
                (
                    Thermo::blendedVolScalarFieldProperty
                    (
                        "cSqr",
                        sqr(dimVelocity),
                        &Thermo::thermoType1::cSqr,
                        &Thermo::thermoType2::cSqr,
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


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<Thermo>::speedOfSound(const label patchi) const
{
    return sqrt
    (
        max
        (
            Thermo::blendedPatchFieldProperty
            (
                &Thermo::thermoType1::cSqr,
                &Thermo::thermoType2::cSqr,
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


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<Thermo>::calcP(const label patchi) const
{
    return
        max
        (
            Thermo::blendedPatchFieldProperty
            (
                &Thermo::thermoType1::p,
                &Thermo::thermoType2::p,
                patchi,
                this->rho_.boundaryField()[patchi],
                this->e_.boundaryField()[patchi],
                this->T_.boundaryField()[patchi]
            ),
            small
        );
}


template<class Thermo>
Foam::scalar Foam::detonatingFluidThermo<Thermo>::calcPi(const label celli) const
{
    const scalar& x = this->xi(celli);
    scalar pi;
   if (x < small)
    {
        pi =
            Thermo::thermoType1::p
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
    }
    else if ((1.0 - x) < small)
    {
        pi =
            Thermo::thermoType2::p
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
    }
    else
    {
        pi =
            Thermo::thermoType2::p
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            )*this->xi(celli)
          + Thermo::thermoType1::p
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            )*(1.0 - this->xi(celli));
    }
    return max(pi, small);
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
    if (this->rho_.time().timeIndex() == 0)
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
