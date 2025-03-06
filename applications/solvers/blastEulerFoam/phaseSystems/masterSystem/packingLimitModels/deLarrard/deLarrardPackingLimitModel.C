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

#include "deLarrardPackingLimitModel.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace packingLimitModels
{
    defineTypeNameAndDebug(deLarrard, 0);

    addToRunTimeSelectionTable
    (
        packingLimitModel,
        deLarrard,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingLimitModels::deLarrard::deLarrard
(
    const dictionary& dict,
    const masterSystem& system
)
:
    binary(dict, system)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingLimitModels::deLarrard::~deLarrard()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::packingLimitModels::deLarrard::aij(const scalar rij) const
{
    return sqrt(1.0 - pow(1.0 - rij, 1.02));
}


Foam::scalar Foam::packingLimitModels::deLarrard::bij(const scalar rij) const
{
    return 1.0 - pow(1.0 - rij, 1.5);
}


// ************************************************************************* //
