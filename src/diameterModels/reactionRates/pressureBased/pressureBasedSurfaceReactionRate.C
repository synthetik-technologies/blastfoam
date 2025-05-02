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

#include "pressureBasedSurfaceReactionRate.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceReactionRates
{
    defineTypeNameAndDebug(pressureBased, 0);
    addToRunTimeSelectionTable(surfaceReactionRate, pressureBased, dictionary);
    addToRunTimeSelectionTable(surfaceReactionRate, pressureBased, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReactionRates::pressureBased::pressureBased(const dictionary& dict)
:
    surfaceReactionRate(dict),
    pScale_(dict.lookup<scalar>("pScale")),
    pExponent_("pExponent", dimless, dict),
    pCoeff_("pCoeff", pow(dimPressure, -pExponent_)*dimLength/dimTime, dict),
    pMin_("pMin", dimPressure, dict.lookupOrDefault<scalar>("pMin", 0.0)),
    offset_("offset", dimLength/dimTime, dict.lookupOrDefault<scalar>("offset", 0.0))
{}


Foam::surfaceReactionRates::pressureBased::pressureBased
(
    const Time& runTime,
    const dictionary& dict
)
:
    pressureBased(dict)
{}


Foam::surfaceReactionRates::pressureBased::pressureBased
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pressureBased(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceReactionRates::pressureBased::~pressureBased()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::surfaceReactionRates::pressureBased::k
(
    const scalar p,
    const scalar T,
    const label
) const
{
    if (p < pMin_.value())
    {
        return 0.0;
    }
    scalar K = pCoeff_.value();
    if (mag(pExponent_.value()) > vSmall)
    {
        K *= pow(p*pScale_, pExponent_.value());
    }
    return offset_.value() + K;
}


Foam::tmp<Foam::volScalarField> Foam::surfaceReactionRates::pressureBased::k
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpk
    (
        volScalarField::New
        (
            typeName + ":k",
            p.mesh(),
            pCoeff_
        )
    );
    volScalarField& K = tmpk.ref();
    if (mag(pExponent_.value()) > vSmall)
    {
        K *= pow(p*pScale_, pExponent_);
    }
    if (offset_.value() > vSmall)
    {
        K += offset_;
    }
    return tmpk*pos(p - pMin_);
}


Foam::tmp<Foam::scalarField> Foam::surfaceReactionRates::pressureBased::k
(
    const fvPatchScalarField& p,
    const fvPatchScalarField& T
) const
{
    tmp<scalarField> tmpk(new scalarField(p.size(), pCoeff_.value()));
    scalarField& K = tmpk.ref();
    if (mag(pExponent_.value()) > vSmall)
    {
        forAll(K, fi)
        {
            K[fi] *= pow(p[fi]*pScale_, pExponent_.value());
        }
    }
    if (mag(offset_.value()) > vSmall)
    {
        forAll(K, fi)
        {
            K[fi] += offset_.value();
        }
    }
    return tmpk;
}

// ************************************************************************* //
