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

#include "YuStandishPackingLimitModel.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingLimitModels::YuStandish::YuStandish
(
    const dictionary& dict,
    const masterSystem& system
)
:
    packingLimitModel(dict, system)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingLimitModels::YuStandish::~YuStandish()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::packingLimitModels::YuStandish::alphaMax
(
    const label celli,
    SortableList<scalar>& ds,
    const bool fixedD
) const
{
    scalar alphap = system_.alpha()[celli];

    if(alphap < system_.residualAlpha().value())
    {
        return minAlphaMax_;
    }

    if (!fixedD)
    {
        ds.reverseSort();
    }

    const UPtrList<phaseModel>& phases(system_.phases());

    scalar maxAlpha = 1.0;

    forAll(ds, i)
    {
        const label phasei = ds.indices()[i];
        const phaseModel& phase1 = phases[phasei];
        scalar alpha1 = phase1[celli];
        if (alpha1 < phase1.residualAlpha().value())
        {
            continue;
        }

        scalar alphaMax1 = phase1.alphaMax();
        scalar d1 = ds[i];

        scalar cxi = alpha1/max(alphap, system_.residualAlpha().value());

        scalar sum = 0.0;

        forAll(ds, j)
        {
            if (i != j)
            {
                const label phasej = ds.indices()[j];
                const phaseModel& phase2 = phases[phasej];
                scalar cxj = phase2[celli]/max(alphap, system_.residualAlpha().value());
                if (phase2[celli] > phase2.residualAlpha().value())
                {
                    continue;
                }
                scalar d2 = ds[j];

                scalar rij = i > j ? d1/d2 : d2/d1;
                if (mag(rij - 1.0) > small)
                {
                    scalar Xij =
                        i > j
                    ? (1.0 - sqr(rij))/(2.0 - alphaMax1)
                    : 1.0 - (1.0 - sqr(rij))/(2.0 - alphaMax1);
                    scalar pij = alphaMax1;

                    if (rij <= 0.741)
                    {
                        pij +=
                            alphaMax1
                        *(1.0 - alphaMax1)
                        *(1.0 - 2.35*rij + 1.35*sqr(rij));
                    }
                    sum += (1.0 - alphaMax1/pij)*cxj/Xij;
                }
            }
        }
        maxAlpha = min(maxAlpha, alphaMax1/(1.0 - sum));
    }
    return max(maxAlpha, minAlphaMax_);
}


// ************************************************************************* //
