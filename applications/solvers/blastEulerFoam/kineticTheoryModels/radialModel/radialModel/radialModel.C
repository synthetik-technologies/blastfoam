/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "radialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(radialModel, 0);

    defineRunTimeSelectionTable(radialModel, dictionary);
}
}

bool Foam::kineticTheoryModels::radialModel::requireKineticTheory() const
{
    if (kt_.valid())
    {
        return true;
    }
    FatalErrorInFunction
        << "Trying to use a kinetic theory based model with a" << nl
        << "non-kinetic theory based system." << endl
        << abort(FatalError);
    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModel::radialModel
(
    const dictionary& dict,
    const masterSystem& system
)
:
    dict_(dict),
    system_(system),
    kt_
    (
        isA<kineticTheorySystem>(system_)
      ? &dynamicCast<const kineticTheorySystem>(system_)
      : nullptr
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModel::~radialModel()
{}


// ************************************************************************* //
