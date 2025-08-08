/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "reactionRateFlameAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reactionRateFlameAreaModel>
Foam::reactionRateFlameAreaModel::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const blastCombustionModel& combModel
)
{
    word reactionRateFlameAreaModelType
    (
        dict.lookup("reactionRateFlameAreaModel")
    );

    Info<< "Selecting reaction rate flame area correlation "
        << reactionRateFlameAreaModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionRateFlameAreaModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown reactionRateFlameAreaModel type "
            << reactionRateFlameAreaModelType << endl << endl
            << "Valid reaction rate flame area types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    const label tempOpen = reactionRateFlameAreaModelType.find('<');

    const word className = reactionRateFlameAreaModelType(0, tempOpen);

    return autoPtr<reactionRateFlameAreaModel>
        (cstrIter()(className, dict, mesh, combModel));
}


// ************************************************************************* //
