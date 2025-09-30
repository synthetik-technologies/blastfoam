/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
------------------------------------------------------------------------*/

#include "surfaceReactionRate.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceReactionRate> Foam::surfaceReactionRate::New
(
    const Time& runTime,
    const dictionary& dict
)
{
    const word reactionRateType(dict.lookup<word>("reactionRate"));

    Info<< "Selecting reactionRate: " << reactionRateType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionRateType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reactionRate type "
            << reactionRateType << endl << endl
            << "Valid reactionRate types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()
    (
        runTime,
        dict.optionalSubDict
        (
            reactionRateType + "ReactionRateCoeffs"
        )
    );
}


Foam::autoPtr<Foam::surfaceReactionRate> Foam::surfaceReactionRate::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word reactionRateType(dict.lookup<word>("reactionRate"));

    Info<< "Selecting fvMesh reactionRate: " << reactionRateType << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(reactionRateType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fvMesh reactionRate type "
            << reactionRateType << endl << endl
            << "Valid fvMesh reactionRate types are : " << endl
            << fvMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()
    (
        mesh,
        dict.optionalSubDict
        (
            reactionRateType + "ReactionRateCoeffs"
        )
    );

}


// ************************************************************************* //
