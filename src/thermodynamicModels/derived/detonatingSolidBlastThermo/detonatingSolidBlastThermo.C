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

#include "detonatingSolidBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastSolidMixture<BasicMixture, ThermoType1, ThermoType2>::HE
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Es(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Es(rho, e, T);
    }
    return thermo2_.Es(rho, e, T)*xi_ + thermo1_.Es(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastSolidMixture<BasicMixture, ThermoType1, ThermoType2>::TRhoE
(
    const scalar T,
    const scalar rho,
    const scalar e
) const
{
    if (xi_ < small)
    {
        return thermo1_.TRhoE(T, rho, e);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.TRhoE(T, rho, e);
    }

    return
        thermo2_.TRhoE(T, rho, e)*xi_
      + thermo1_.TRhoE(T, rho, e)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastSolidMixture<BasicMixture, ThermoType1, ThermoType2>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Cp(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Cp(rho, e, T);
    }
    return thermo2_.Cp(rho, e, T)*xi_ + thermo1_.Cp(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastSolidMixture<BasicMixture, ThermoType1, ThermoType2>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.Cv(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.Cv(rho, e, T);
    }
    return thermo2_.Cv(rho, e, T)*xi_ + thermo1_.Cv(rho, e, T)*(1.0 - xi_);
}


template<class BasicMixture, class ThermoType1, class ThermoType2>
Foam::scalar
Foam::detonatingBlastSolidMixture<BasicMixture, ThermoType1, ThermoType2>::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (xi_ < small)
    {
        return thermo1_.kappa(rho, e, T);
    }
    else if ((1.0 - xi_) < small)
    {
        return thermo2_.kappa(rho, e, T);
    }
    return thermo2_.kappa(rho, e, T)*xi_ + thermo1_.kappa(rho, e, T)*(1.0 - xi_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingSolidBlastThermo<Thermo>::detonatingSolidBlastThermo
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
        dict.subDict("reactants"),
        dict.subDict("products"),
        phaseName,
        masterName
    ),
    activation_
    (
        activationModel::New
        (
            mesh,
            dict,
            phaseName
        )
    ),
    afterburn_
    (
        afterburnModel::New
        (
            mesh,
            dict,
            phaseName
        )
    ),
    mixture_(*this, *this, activation_())
{
    this->initializeFields();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::initializeModels()
{
    activation_->initializeModels();
    afterburn_->initializeModels();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingSolidBlastThermo<Thermo>::~detonatingSolidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::update()
{
    activation_->update();
    afterburn_->update();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::solve()
{
    activation_->solve();
    afterburn_->solve();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::postUpdate()
{
    activation_->postUpdate();
    afterburn_->postUpdate();
}



template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::updateRho()
{
    volScalarField& rhoRef(this->rho_);
    volScalarField rhoNew
    (
        Thermo::blendedVolScalarFieldProperty
        (
            "rho",
            dimDensity,
            &Thermo::thermoType1::rho0,
            &Thermo::thermoType2::rho0
        )
    );
    forAll(rhoRef, celli)
    {
        rhoRef[celli] = rhoNew[celli];
    }
    forAll(rhoRef.boundaryField(), patchi)
    {
        forAll(rhoRef.boundaryField()[patchi], facei)
        {
            rhoRef.boundaryFieldRef()[patchi][facei] =
                rhoNew.boundaryField()[patchi][facei];
        }
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingSolidBlastThermo<Thermo>::calce() const
{
    tmp<volScalarField> e
    (
        Thermo::blendedVolScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &Thermo::thermoType1::Es,
            &Thermo::thermoType2::Es,
            this->rho_,
            this->e_,
            this->T_
        )
    );

    //- Add detonation energy to initially reacted material
    if (this->rho_.time().timeIndex() == 0)
    {
        e.ref() += activation_->e0()*activation_->lambda();
    }

    return e;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingSolidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (activation_->ESource() + afterburn_->ESource())*this->rho_
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::detonatingSolidBlastThermo<Thermo>::kappa() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &Thermo::thermoType1::kappa,
        &Thermo::thermoType2::kappa,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volVectorField>
Foam::detonatingSolidBlastThermo<Thermo>::Kappa() const
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
    const scalarField& rhoCells = this->rho_;
    const scalarField& eCells = this->e_;
    const scalarField& TCells = this->T_;

    forAll(KappaCells, celli)
    {
        scalar x = cellx(celli);
        if (x < small)
        {
            Kappa[celli] =
                Thermo::thermoType1::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                );
        }
        else if ((1.0 - x) < small)
        {
            Kappa[celli] =
                Thermo::thermoType2::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                );
        }
        else
        {
            Kappa[celli] =
                Thermo::thermoType2::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                )*x
              + Thermo::thermoType1::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                )*(1.0 - x);
        }
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        vectorField& Kappap = KappaBf[patchi];
        const scalarField& pRho = this->rho_.boundaryField()[patchi];
        const scalarField& pe = this->e_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        tmp<scalarField> xp(this->x(patchi));

        forAll(Kappap, facei)
        {
            const scalar& x = xp()[facei];
            if (x < small)
            {
                Kappap[facei] =
                    Thermo::thermoType1::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    );
            }
            else if ((1.0 - x) < small)
            {
                Kappap[facei] =
                    Thermo::thermoType2::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    );
            }
            else
            {
                Kappap[facei] =
                    Thermo::thermoType2::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    )*x
                  + Thermo::thermoType1::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    )*(1.0 - x);
            }
        }
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::vectorField>
Foam::detonatingSolidBlastThermo<Thermo>::Kappa(const label patchi) const
{
    const scalarField& pRho = this->rho_.boundaryField()[patchi];
    const scalarField& pe = this->e_.boundaryField()[patchi];
    const scalarField& pT = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pe.size()));

    vectorField& Kappap = tKappa.ref();
    tmp<scalarField> xp(this->x(patchi));

    forAll(pe, facei)
    {
        const scalar& x = xp()[facei];
        if (x < small)
        {
            Kappap[facei] =
                Thermo::thermoType1::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                );
        }
        else if ((1.0 - x) < small)
        {
            Kappap[facei] =
                Thermo::thermoType2::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                );
        }
        else
        {
            Kappap[facei] =
                Thermo::thermoType2::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                )*x
                + Thermo::thermoType1::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                )*(1.0 - x);
        }
    }

    return tKappa;
}


// ************************************************************************* //
