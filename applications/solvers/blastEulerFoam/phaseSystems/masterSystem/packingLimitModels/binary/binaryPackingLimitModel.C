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

#include "binaryPackingLimitModel.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(binary, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingLimitModels::binary::binary
(
    const dictionary& dict,
    const masterSystem& system
)
:
    packingLimitModel(dict, system)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingLimitModels::binary::~binary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::packingLimitModels::binary::alphaMax
(
    const label celli,
    SortableList<scalar>& ds,
    const bool fixedD
) const
{
    const scalar alphap = system_.alpha()[celli];
    if(alphap < system_.residualAlpha().value())
    {
        return minAlphaMax_;
    }

    if (!fixedD)
    {
        ds.reverseSort();
    }

    const UPtrList<phaseModel>& phases(system_.phases());

    label nValid = 0;
    scalar maxAlpha = 1.0;
    forAll(ds, i)
    {
        const label phasei = ds.indices()[i];
        const phaseModel& phase1 = phases[phasei];
        const scalar alpha1 = phase1[celli];
        if (alpha1 < phase1.residualAlpha().value())
        {
            continue;
        }
        nValid++;

        const scalar alphaMax1 = phase1.alphaMax();
        const scalar d1 = ds[i];

        scalar denom = 0.0;

        for (label j = 0; j < i; j++)
        {
            const label phasej = ds.indices()[j];
            const phaseModel& phase2 = phases[phasej];
            const scalar alpha2 = phase2[celli];
            if (alpha2 < phase2.residualAlpha().value())
            {
                continue;
            }
            const scalar alphaMax2 = phase2.alphaMax();

            const scalar rij = d1/ds[j];
            denom +=
                (
                    1.0
                  - alphaMax1
                  + this->bij(rij)*alphaMax1*(1.0 - 1.0/alphaMax2)
                )*alpha2/alphap;
        }
        for (label j = i+1; j < ds.size(); j++)
        {
            const label phasej = ds.indices()[j];
            const phaseModel& phase2 = phases[phasej];
            const scalar alpha2 = phase2[celli];
            if (alpha2 < phase2.residualAlpha().value())
            {
                continue;
            }
            const scalar alphaMax2 = phase2.alphaMax();

            const scalar rij = ds[j]/d1;
            denom +=
                (
                    1.0
                  - this->aij(rij)*alphaMax1/alphaMax2
                )*alpha2/alphap;

        }
        maxAlpha = min(maxAlpha, alphaMax1/(1.0 - denom));
    }

    if (!nValid)
    {
        return minAlphaMax_;
    }
    return max(maxAlpha, minAlphaMax_);
}


// ************************************************************************* //
