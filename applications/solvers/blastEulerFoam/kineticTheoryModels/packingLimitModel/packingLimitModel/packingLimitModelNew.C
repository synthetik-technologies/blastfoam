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

#include "packingLimitModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::kineticTheoryModels::packingLimitModel>
Foam::kineticTheoryModels::packingLimitModel::New
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
{
    word packingLimitModelType
    (
        dict.lookupOrDefault<word>("packingLimitModel", "constant")
    );

    Info<< "Selecting packingLimitModel "
        << packingLimitModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(packingLimitModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "packingLimitModel::New(const dictionary&) : " << endl
            << "    unknown packingLimitModelType type "
            << packingLimitModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid packingLimitModelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<packingLimitModel>(cstrIter()(dict, kt));
}


// ************************************************************************* //
