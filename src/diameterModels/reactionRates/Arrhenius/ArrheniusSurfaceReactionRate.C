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

#include "ArrheniusSurfaceReactionRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceReactionRates
{
    defineTypeNameAndDebug(Arrhenius, 0);
    addToRunTimeSelectionTable(surfaceReactionRate, Arrhenius, dictionary);
    addToRunTimeSelectionTable(surfaceReactionRate, Arrhenius, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReactionRates::Arrhenius::Arrhenius(const dictionary& dict)
:
    surfaceReactionRate(dict),
    A_("A", inv(dimTime), dict),
    beta_("beta", dimless, dict),
    Ta_("Ta", dimTemperature, dict)
{}


Foam::surfaceReactionRates::Arrhenius::Arrhenius
(
    const Time& runTime,
    const dictionary& dict
)
:
    Arrhenius(dict)
{}


Foam::surfaceReactionRates::Arrhenius::Arrhenius
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    Arrhenius(dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceReactionRates::Arrhenius::~Arrhenius()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::surfaceReactionRates::Arrhenius::k
(
    const scalar p,
    const scalar T,
    const label
) const
{
    scalar K = A_.value();

    if (mag(beta_.value()) > vSmall)
    {
        K *= pow(T, beta_.value());
    }
    if (mag(Ta_.value()) > vSmall)
    {
        K *= exp(-Ta_.value()/T);
    }
    return K;
}


Foam::tmp<Foam::volScalarField> Foam::surfaceReactionRates::Arrhenius::k
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
                typeName + ":k",
                p.time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p.mesh(),
            A_
        )
    );
    volScalarField& K = tmpk.ref();
    if (mag(beta_.value()) > vSmall)
    {
        K *= pow(T, beta_);
    }
    if (mag(Ta_.value()) > vSmall)
    {
        K *= exp(-Ta_/T);
    }
    return tmpk;
}


Foam::tmp<Foam::scalarField> Foam::surfaceReactionRates::Arrhenius::k
(
    const fvPatchScalarField& p,
    const fvPatchScalarField& T
) const
{
    tmp<scalarField> tmpk(new scalarField(p.size(), A_.value()));
    scalarField& K = tmpk.ref();
    if (mag(beta_.value()) > vSmall)
    {
        forAll(K, fi)
        {
            K[fi] *= pow(T[fi], beta_.value());
        }
    }
    if (mag(Ta_.value()) > vSmall)
    {
        forAll(K, fi)
        {
            K[fi] *= exp(-Ta_.value()/T[fi]);
        }
    }
    return tmpk;
}


// ************************************************************************* //
