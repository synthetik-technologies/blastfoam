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

#include "tabulatedPReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{
    defineTypeNameAndDebug(tabulatedP, 0);
    addToRunTimeSelectionTable(reactionRate, tabulatedP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRates::tabulatedP::tabulatedP(const dictionary& dict)
:
    reactionRate(dict),
    table_
    (
        dict.found("file")
      ? lookupTable1D
        (
            dict.lookup<fileName>("file"),
            dict.lookupOrDefault<word>("mod", "none"),
            dict.lookupOrDefault<word>("xMod", "none"),
            dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
        )
      : lookupTable1D
        (
            dict.lookup<scalarField>("p"),
            dict.lookup<scalarField>("dxdt"),
            dict.lookupOrDefault<word>("mod", "none"),
            dict.lookupOrDefault<word>("xMod", "none"),
            dict.lookupOrDefault<word>("interpolationScheme", "linearClamp"),
            dict.lookupOrDefault("correct", false)
        )
    ),
    pScale_(dict.lookup<scalar>("pScale"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRates::tabulatedP::~tabulatedP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::reactionRates::tabulatedP::k
(
    const scalar& p,
    const scalar& T
) const
{
    return table_.lookup(p*pScale_);
}


Foam::tmp<Foam::volScalarField> Foam::reactionRates::tabulatedP::k
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpk
    (
        new volScalarField
        (
            IOobject
            (
                "tabulatedP:k",
                p.time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p.mesh(),
            dimensionedScalar(dimLength/dimTime, 0.0)
        )
    );
    volScalarField& K = tmpk.ref();
    forAll(K, celli)
    {
        K[celli] = table_.lookup(p[celli]*pScale_);
    }
    return tmpk;
}

// ************************************************************************* //
