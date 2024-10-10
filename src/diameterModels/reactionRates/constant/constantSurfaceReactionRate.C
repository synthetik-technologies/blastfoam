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

#include "constantSurfaceReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceReactionRates
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(surfaceReactionRate, constant, dictionary);
    addToRunTimeSelectionTable(surfaceReactionRate, constant, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReactionRates::constant::constant(const dictionary& dict)
:
    surfaceReactionRate(dict),
    rate_("rate", inv(dimTime), dict)
{}


Foam::surfaceReactionRates::constant::constant
(
    const Time& runTime,
    const dictionary& dict
)
:
    constant(dict)
{}


Foam::surfaceReactionRates::constant::constant
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    constant(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceReactionRates::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::surfaceReactionRates::constant::k
(
    const scalar p,
    const scalar T,
    const label
) const
{
    return rate_.value();
}


Foam::tmp<Foam::volScalarField> Foam::surfaceReactionRates::constant::k
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarField::New
    (
        typeName + ":k",
        p.mesh(),
        rate_
    );
}

// ************************************************************************* //
