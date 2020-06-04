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

#include "basicSolidThermoModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicSolidThermoModel<Thermo>::basicSolidThermoModel
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master
)
:
    Thermo
    (
        name,
        mesh,
        dict,
        master
    )
{
    volScalarField& rhoRef(solidThermoModel::rho_);
    forAll(rhoRef, celli)
    {
        rhoRef[celli] = Thermo::thermoType::rho0();
    }
    forAll(rhoRef.boundaryField(), patchi)
    {
        forAll(rhoRef.boundaryField()[patchi], facei)
        {
            rhoRef.boundaryFieldRef()[patchi][facei] =
                Thermo::thermoType::rho0();
        }
    }
    this->e_ = this->E();
    this->T_ = this->calcT();
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicSolidThermoModel<Thermo>::~basicSolidThermoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::basicSolidThermoModel<Thermo>::correct()
{
    if (this->master_)
    {
        this->T_ = this->calcT();
    }

    this->Thermo::alpha_ = Thermo::volScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &Thermo::thermoType::kappa,
        this->rho_,
        this->e_,
        this->T_
    )/this->Cv();
}


template<class Thermo>
Foam::tmp<Foam::volVectorField>
Foam::basicSolidThermoModel<Thermo>::Kappa() const
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
    const scalarField& rhoCells = solidThermoModel::rho_;
    const scalarField& eCells = this->e_;
    const scalarField& TCells = this->T_;

    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            Thermo::thermoType::Kappa
            (
                rhoCells[celli],
                eCells[celli],
                TCells[celli]
            );
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        vectorField& Kappap = KappaBf[patchi];
        const scalarField& pRho = solidThermoModel::rho_.boundaryField()[patchi];
        const scalarField& pe = this->e_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] =
                Thermo::thermoType::Kappa(pRho[facei], pe[facei], pT[facei]);
        }
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::vectorField>
Foam::basicSolidThermoModel<Thermo>::Kappa(const label patchi) const
{
    const scalarField& pRho = solidThermoModel::rho_.boundaryField()[patchi];
    const scalarField& pe = this->e_.boundaryField()[patchi];
    const scalarField& pT = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pe.size()));

    vectorField& Kappap = tKappa.ref();

    forAll(pe, facei)
    {
        Kappap[facei] =
            Thermo::thermoType::Kappa(pRho[facei], pe[facei], pT[facei]);
    }

    return tKappa;
}

// ************************************************************************* //
