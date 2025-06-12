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

#include "diameterReactionRate.H"
#include "lengthDiameterReactionRate.H"
#include "surfaceReactionRate.H"
#include "diameterModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterReactionRate> Foam::diameterReactionRate::New
(
    const diameterModel& dModel,
    const dictionary& dict
)
{
    const word reactionRateType(dict.lookup<word>("reactionRate"));

    if
    (
        surfaceReactionRate::dictionaryConstructorTablePtr_->found
        (
            reactionRateType
        )
    )
    {
        return autoPtr<diameterReactionRate>
        (
            new diameterReactionRates::length
            (
                dModel,
                surfaceReactionRate::New(dModel.d().mesh(), dict)
            )
        );
    }
    wordHashSet types(dictionaryConstructorTablePtr_->toc());

    if (dictionaryConstructorTablePtr_)
    {
        Info<< "Selecting reactionRate: " << reactionRateType << endl;
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(reactionRateType);

        if (cstrIter != dictionaryConstructorTablePtr_->end())
        {
            return cstrIter()
            (
                dModel,
                dict.optionalSubDict
                (
                    reactionRateType + "ReactionRateCoeffs"
                )
            );
        }
        else
        {
            types += surfaceReactionRate::dictionaryConstructorTablePtr_->toc();
        }
    }
    FatalErrorInFunction
        << "Unknown reactionRate type "
        << reactionRateType << endl << endl
        << "Valid reactionRate types are : " << endl
        << types.toc()
        << exit(FatalError);

    return autoPtr<diameterReactionRate>();
}


// ************************************************************************* //
