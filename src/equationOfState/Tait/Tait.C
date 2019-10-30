/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "Tait.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace equationOfStates
{
    defineTypeNameAndDebug(Tait, 0);
    addToRunTimeSelectionTable(equationOfState, Tait, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationOfStates::Tait::Tait
(
    const volScalarField& e,
    const dictionary& dict
)
:
    equationOfState(e, dict),
    a_("a", dimensionSet(1, -1, -2, 0, 0, 0, 0), dict),
    b_("b", dimensionSet(1, -1, -2, 0, 0, 0, 0), dict),
    gamma_("gamma", dimensionSet(0, 0, 0, 0, 0, 0, 0), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationOfStates::Tait::~Tait()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::Tait::Gamma() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Gamma",
                rho_.time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho_.mesh(),
            gamma_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::Tait::Pi() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Pi",
                rho_.time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho_.mesh(),
            gamma_*(b_ - a_)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::equationOfStates::Tait::delta() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "delta",
                rho_.time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho_.mesh(),
            dimensionedScalar("delta", sqr(dimVelocity), 0.0)
        )
    );
}


// ************************************************************************* //
