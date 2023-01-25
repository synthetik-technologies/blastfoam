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
namespace kineticTheoryModels
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::constant::constant
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    packingLimitModel(dict, kt),
    maxAlpha_
    (
        dict.lookupOrDefault
        (
            "alphaMax",
            -1
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::packingLimitModels::constant::alphaMax() const
{
    const UPtrList<phaseModel>& phases(kt_.phases());

    scalar alphaMax = maxAlpha_;
    if (alphaMax < 0)
    {
        forAll(phases, phasei)
        {
            alphaMax = max(phases[phasei].alphaMax(), alphaMax);
        }
    }
    return
        volScalarField::New
        (
            "alphaMax",
            mesh_,
            alphaMax,
            zeroGradientFvPatchScalarField::typeName
        );
}

Foam::scalar
Foam::kineticTheoryModels::packingLimitModels::constant::alphaMax
(
    const label celli,
    const SortableList<scalar>& ds
) const
{
    NotImplemented;
    return maxAlpha_;
}


bool Foam::kineticTheoryModels::packingLimitModels::constant::read()
{
    maxAlpha_ = dict_.lookup<scalar>("alphaMax");

    return true;
}


// ************************************************************************* //
