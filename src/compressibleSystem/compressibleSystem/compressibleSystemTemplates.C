/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "compressibleSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class SelectionTable>
Foam::autoPtr<Foam::compressibleSystem> Foam::compressibleSystem::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& defaultType,
    SelectionTable* tablePtr
)
{
    typename SelectionTable::iterator cstrIter = tablePtr->begin();
    if (tablePtr->size() > 1)
    {
        word compressibleSystemType
        (
            dict.lookupOrDefault<word>("systemType", defaultType)
        );

        Info<< "Selecting " << compressibleSystemType
            << " compressibleSystem" << endl;

        cstrIter = tablePtr->find(compressibleSystemType);

        if (cstrIter == tablePtr->end())
        {
            FatalErrorInFunction
                << "Unknown compressibleSystem type "
                << compressibleSystemType << endl << endl
                << "Valid compressibleSystem types are : " << endl
                << tablePtr->sortedToc()
                << exit(FatalError);
        }
    }

    return cstrIter()(mesh);
}


// ************************************************************************* //
