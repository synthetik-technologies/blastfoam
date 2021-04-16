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

#include "pressureBasedReactionRate.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{
    defineTypeNameAndDebug(pressureBased, 0);
    addToRunTimeSelectionTable(reactionRate, pressureBased, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRates::pressureBased::pressureBased(const dictionary& dict)
:
    reactionRate(dict),
    pScale_(dict.lookup<scalar>("pScale")),
    pExponent_("pExponent", dimless, dict),
    pCoeff_("pCoeff", pow(dimPressure, -pExponent_)*dimLength/dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRates::pressureBased::~pressureBased()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::reactionRates::pressureBased::k
(
    const scalar& p,
    const scalar& T
) const
{
    scalar K = pCoeff_.value();
    if (mag(pExponent_.value()) > vSmall)
    {
        K *= pow(p*pScale_, pExponent_.value());
    }
    return K;
}


Foam::tmp<Foam::volScalarField> Foam::reactionRates::pressureBased::k
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
                "pressureBased:k",
                p.time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p.mesh(),
            pCoeff_
        )
    );
    volScalarField& K = tmpk.ref();
    if (mag(pExponent_.value()) > vSmall)
    {
        K *= pow(p*pScale_, pExponent_);
    }
    return tmpk;
}

// ************************************************************************* //
