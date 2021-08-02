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

#include "detonatingSolidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingSolidThermo<Thermo>::detonatingSolidThermo
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    Thermo
    (
        phaseName,
        mesh,
        dict,
        dict.subDict("reactants"),
        dict.subDict("products"),
        master,
        masterName
    )
{}


template<class Thermo>
void Foam::detonatingSolidThermo<Thermo>::initializeModels()
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

    updateRho();

    this->e_ = this->E();

    //- Add detonation energy
    if
    (
        this->rho_.time().value()
     == this->rho_.time().startTime().value()
    )
    {
        this->e_ += activation_->e0()*activation_->lambda();
    }

    this->T_ = this->calcT();
    this->correct();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingSolidThermo<Thermo>::~detonatingSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingSolidThermo<Thermo>::update()
{
    activation_->update();
    afterburn_->update();
}


template<class Thermo>
void Foam::detonatingSolidThermo<Thermo>::solve()
{
    activation_->solve();
    afterburn_->solve();
}


template<class Thermo>
void Foam::detonatingSolidThermo<Thermo>::postUpdate()
{
    activation_->postUpdate();
    afterburn_->postUpdate();
}



template<class Thermo>
void Foam::detonatingSolidThermo<Thermo>::updateRho()
{
    volScalarField& rhoRef(solidBlastThermo::rho_);
    volScalarField rhoNew
    (
        Thermo::blendedVolScalarFieldProperty
        (
            IOobject::groupName("rhoNew", basicBlastThermo::name_),
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
void Foam::detonatingSolidThermo<Thermo>::correct()
{
    if (this->master_)
    {
        this->T_ = this->calcT();
    }

    this->Thermo::alpha_ = Thermo::blendedVolScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &Thermo::thermoType1::kappa,
        &Thermo::thermoType2::kappa,
        this->rho_,
        this->e_,
        this->T_
    )/this->Cv();

    // Update density
    volScalarField& rhoRef(solidBlastThermo::rho_);
    rhoRef =
        Thermo::blendedVolScalarFieldProperty
        (
            IOobject::groupName("rhos", basicBlastThermo::name_),
            dimDensity,
            &Thermo::thermoType1::rho0,
            &Thermo::thermoType2::rho0
        );
    rhoRef.correctBoundaryConditions();
}

template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingSolidThermo<Thermo>::ESource() const
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

template<class Thermo>
Foam::tmp<Foam::volVectorField>
Foam::detonatingSolidThermo<Thermo>::Kappa() const
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
    const scalarField& rhoCells = solidBlastThermo::rho_;
    const scalarField& eCells = this->e_;
    const scalarField& TCells = this->T_;

    forAll(KappaCells, celli)
    {
        scalar x = xi(celli);
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
        const scalarField& pRho = solidBlastThermo::rho_.boundaryField()[patchi];
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
Foam::detonatingSolidThermo<Thermo>::Kappa(const label patchi) const
{
    const scalarField& pRho = solidBlastThermo::rho_.boundaryField()[patchi];
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
