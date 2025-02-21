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

#include "constantPackingLimit.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(constant, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        constant,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingLimitModels::constant::constant
(
    const dictionary& dict,
    const masterSystem& system
)
:
    packingLimitModel(dict, system),
    maxAlpha_
    (
        dict.lookupOrDefault
        (
            "alphaMax",
            -1.0
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingLimitModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::packingLimitModels::constant::updateAlphaMax
(
    volScalarField::Internal& alphaMaxI
) const
{
    const UPtrList<phaseModel>& phases(system_.phases());

    if (maxAlpha_ < 0)
    {
        forAll(phases, phasei)
        {
            maxAlpha_ = max(phases[phasei].alphaMax(), maxAlpha_);
        }
    }
    alphaMaxI = maxAlpha_;
}

Foam::scalar Foam::packingLimitModels::constant::alphaMax
(
    const label celli,
    const SortableList<scalar>& ds
) const
{
    return maxAlpha_;
}


// ************************************************************************* //
