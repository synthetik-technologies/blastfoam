/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "YuStandishPackingLimit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(YuStandish, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        YuStandish,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::YuStandish::YuStandish
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    packingLimitModel(dict, kt),
    residualAlpha_(dict_.lookupType<scalar>("residualAlpha"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::YuStandish::~YuStandish()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::kineticTheoryModels::packingLimitModels::YuStandish::alphaMax
(
    const label celli,
    const scalarList& ds
) const
{
    scalar alphap = kt_.alphap()[celli];

    if(alphap < kt_.residualAlpha().value())
    {
        return kt_.minAlphaMax();
    }

    scalar maxAlpha = 1.0;

    forAll(ds, phasei)
    {
        const phaseModel& phase1 = kt_.fluid().phases()[phasei];
        scalar alpha1 = phase1[celli];

        scalar alphaMax1 = phase1.alphaMax();
        scalar d1 = ds[phasei];

        scalar cxi = alpha1/max(alphap, residualAlpha_);

        scalar sum = 0.0;

        forAll(ds, phasej)
        {
            if (phasej != phasei)
            {
                scalar d2 = ds[phasej];

                scalar rij = d1/d2;
                scalar Xij = 1.0;
                scalar pij = alphaMax1;

                if (rij >= 1)
                {
                    rij = 1.0/rij;
                    Xij = (1.0 - sqr(rij))/(2.0 - alphaMax1);
                }
                else
                {
                    rij = d2/d1;
                    Xij = 1.0 - (1.0 - sqr(rij))/(2.0 - alphaMax1);
                }

                if (rij <= 0.741)
                {
                    pij +=
                        alphaMax1
                       *(1.0 - alphaMax1)
                       *(1.0 - 2.35*rij + 1.35*sqr(rij));
                }
                sum += (1.0 - alphaMax1/pij)*cxi/Xij;
            }
        }
        maxAlpha = min(maxAlpha, alphaMax1/(1.0 - sum));
    }

    return maxAlpha;
}


bool Foam::kineticTheoryModels::packingLimitModels::YuStandish::read()
{
    return true;
}


// ************************************************************************* //
