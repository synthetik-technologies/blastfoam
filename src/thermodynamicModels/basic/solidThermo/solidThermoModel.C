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

#include "solidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(solidThermoModel, 0);
    defineRunTimeSelectionTable(solidThermoModel, basicSolid);
    defineRunTimeSelectionTable(solidThermoModel, detonatingSolid);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermoModel::solidThermoModel
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    basicThermoModel
    (
        name,
        mesh,
        dict,
        master,
        masterName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermoModel::~solidThermoModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidThermoModel::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            e_.mesh(),
            dimensionedScalar("0", dimDynamicViscosity, 0.0)
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::solidThermoModel::mu(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField(this->rho_.boundaryField()[patchi].size(), 0.0)
    );
}


Foam::tmp<Foam::volScalarField> Foam::solidThermoModel::nu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->rho_.mesh(),
            dimensionedScalar("0", dimViscosity, 0.0)
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::solidThermoModel::nu(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField(this->rho_.boundaryField()[patchi].size(), 0.0)
    );
}


Foam::scalar Foam::solidThermoModel::nui(const label) const
{
    return 0.0;
}

// ************************************************************************* //
