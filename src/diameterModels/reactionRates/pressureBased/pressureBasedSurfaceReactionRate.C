/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
     \\/     M anipulation  | Synthetik Applied Technologies
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
    pExponent_("pExponent", dimless, dict),
    pCoeff_
    (
        "pCoeff",
        pow(dimPressure, -pExponent_)*dimLength/dimTime,
        dict.lookup<scalar>("pCoeff") // No units, handled later
    ),
    pMin_
    (
        "pMin",
        dimPressure,
        dict.lookupOrDefault<scalar>("pMin", dimPressure, 0.0)
    ),
    offset_
    (
        "offset",
        dimLength/dimTime,
        dict.lookupOrDefault<scalar>("offset", dimLength, 0.0)
    )
{
    // Convert pressure coefficient
    pCoeff_ *= setPCoeffUnits(dict, pExponent_.value());
}


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
        K *= pow(p, pExponent_.value());
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
        K *= pow(p, pExponent_);
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
        K *= pow(p, pExponent_.value());
    }
    if (mag(offset_.value()) > vSmall)
    {
        K += offset_.value();
    }
    if (mag(pMin_.value()) > vSmall)
    {
        K *= pos(p - pMin_.value());
    }

    return tmpk;
}

Foam::scalar
Foam::surfaceReactionRates::pressureBased::setPCoeffUnits
(
    const dictionary& dict,
    const scalar pExponent
)
{
    bool foundPScale = dict.found("pScale");

    if (foundPScale && dict.found("pCoeffUnits"))
    {
        FatalErrorInFunction
            << "Both \"pScale\" and \"pCoeffUnits\" were specified."
            << endl
            << "User must only specify one for proper conversion!"
            << endl
            << abort(FatalError);
    }

    if (foundPScale)
    {
        scalar pScale = dict.lookup<scalar>("pScale");
        return pow(pScale, pExponent);
    }
    else
    {
        // Convert units
        PtrList<unitConversion> pCoeffUnits;

        ITstream& is = dict.lookup("pCoeffUnits");
        while (is.good())
        {
            // Construct directly from stream
            pCoeffUnits.append(new unitConversion(is));
        }

        scalar factor = 1.0;

        if (pCoeffUnits.size() > 2)
        {
            FatalIOErrorInFunction(dict)
                << "Only pressure and velocity unit conversions can be provided, "
                << "but found "
                << pCoeffUnits << endl
                << exit(FatalIOError);
        }
        else if (pCoeffUnits.size())
        {
            bool setPressure = false;
            bool setVelocity = false;
            forAll(pCoeffUnits, j)
            {
                const unitConversion& conv = pCoeffUnits[j];

                if (conv.dimensions() == dimPressure && !setPressure)
                {
                    factor *= pow(conv.toStandard(1.0), -pExponent);
                    setPressure = true;
                }
                else if (conv.dimensions() == dimVelocity && !setVelocity)
                {
                    factor *= conv.toStandard(1.0);
                    setVelocity = true;
                }
                else
                {
                    FatalIOErrorInFunction(dict)
                        << "Only pressure or velocity units can be used, but found "
                        << conv.dimensions() << endl
                        << exit(FatalIOError);
                }

                if (setVelocity && setPressure) break;
            }
        }

        return factor;
    }
}

// ************************************************************************* //
