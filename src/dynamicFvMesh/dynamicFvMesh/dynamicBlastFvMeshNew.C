/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "staticBlastFvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicBlastFvMesh> Foam::dynamicBlastFvMesh::New(const IOobject& io)
{
    IOobject dictHeader(dynamicMeshDictIOobject(io));

    if (dictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictHeader);

        const word dynamicFvMeshTypeName
        (
            dict.lookup<word>("dynamicFvMesh")
        );

        Info<< "Selecting dynamicBlastFvMesh " << dynamicFvMeshTypeName << endl;

        libs.open
        (
            dict,
            "blastDynamicFvMeshLibs",
            dictionaryConstructorTablePtr_
        );

        if (!dictionaryConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "dynamicBlastFvMesh table is empty"
                << exit(FatalError);
        }

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(dynamicFvMeshTypeName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown dynamicBlastFvMesh type "
                << dynamicFvMeshTypeName << nl << nl
                << "Valid dynamicBlastFvMesh types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<dynamicBlastFvMesh>(cstrIter()(io));
    }
    else
    {
        return autoPtr<dynamicBlastFvMesh>(new staticBlastFvMesh(io));
    }
}


// ************************************************************************* //
