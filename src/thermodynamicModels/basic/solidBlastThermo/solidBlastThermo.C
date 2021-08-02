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

#include "solidBlastThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(solidBlastThermo, 0);
    defineRunTimeSelectionTable(solidBlastThermo, basicSolid);
    defineRunTimeSelectionTable(solidBlastThermo, detonatingSolid);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBlastThermo::solidBlastThermo
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    basicBlastThermo
    (
        name,
        mesh,
        dict,
        master,
        masterName
    )
{}


Foam::solidBlastThermo::solidBlastThermo
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    basicBlastThermo
    (
        phaseName,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    )
{}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBlastThermo::~solidBlastThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidBlastThermo::mu() const
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


Foam::tmp<Foam::scalarField> Foam::solidBlastThermo::mu(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField(this->rho_.boundaryField()[patchi].size(), 0.0)
    );
}


Foam::tmp<Foam::volScalarField> Foam::solidBlastThermo::nu() const
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


Foam::tmp<Foam::scalarField> Foam::solidBlastThermo::nu(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField(this->rho_.boundaryField()[patchi].size(), 0.0)
    );
}


Foam::scalar Foam::solidBlastThermo::nui(const label) const
{
    return 0.0;
}

// ************************************************************************* //
