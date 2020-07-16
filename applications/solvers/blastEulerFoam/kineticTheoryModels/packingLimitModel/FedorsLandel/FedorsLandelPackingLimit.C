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

#include "FedorsLandelPackingLimit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(FedorsLandel, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        FedorsLandel,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::FedorsLandel::FedorsLandel
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    packingLimitModel(dict, kt),
    residualAlpha_(dict_.lookupType<scalar>("residualAlpha"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::FedorsLandel::~FedorsLandel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::kineticTheoryModels::packingLimitModels::FedorsLandel::alphaMax
(
    const label celli,
    const scalarList& ds
) const
{
    if (ds.size() != 2)
    {
        FatalErrorInFunction
            << typeName << " packing limit model only supports bi-disperse "
            << "particle distributions, " << ds.size() << " in use."
            << exit(FatalError);
    }

     scalar alphap = kt_.alphap()[celli];

    if(alphap < kt_.fluid().phases()[0].residualAlpha().value())
    {
        return kt_.fluid().phases()[0].alphaMax();
    }


    const phaseModel& phase1 = kt_.fluid().phases()[0];
    scalar alpha1 = phase1[celli];
    scalar alphaMax1 = phase1.alphaMax();
    scalar d1 = ds[0];
    scalar cx1 = alpha1/max(alphap, residualAlpha_);

    scalar alphaMax2 = kt_.fluid().phases()[1].alphaMax();
    scalar d2 = ds[1];
    scalar cx2 = 1.0 - cx1;

    scalar cxMax = alphaMax1/(alphaMax1 + (1.0 - alphaMax1)*alphaMax2);

    scalar r21;
    if (d1 < d2)
    {
        r21 = d1/d2;
    }
    else
    {
        r21 = d2/d1;
    }

    if (cx1 <= cxMax)
    {
        return
        (
            (alphaMax1 - alphaMax2)
          + (1.0 - sqrt(r21))*(1.0 - alphaMax1)*alphaMax2
           *(alphaMax1 + (1.0 - alphaMax1)*alphaMax2)*cx1/alphaMax1
          + alphaMax2
        );
    }
    else
    {
        return
        (
            (1.0 - sqrt(r21))
           *(alphaMax1 + (1.0 - alphaMax1)*alphaMax2)*cx2
          + alphaMax2
        );
    }
}


bool Foam::kineticTheoryModels::packingLimitModels::FedorsLandel::read()
{
    return true;
}


// ************************************************************************* //
