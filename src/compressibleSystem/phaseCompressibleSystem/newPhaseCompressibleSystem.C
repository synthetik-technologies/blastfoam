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

#include "phaseCompressibleSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseCompressibleSystem> Foam::phaseCompressibleSystem::New
(
    const fvMesh& mesh
)
{
    word phaseCompressibleSystemType(word::null);

    // Create temporary phase properties to lookup type
    // not store in the database to remove possible conflict
    Info<< "Reading phaseProperties dictionary\n" << endl;
    IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (phaseProperties.found("systemType"))
    {
        phaseCompressibleSystemType =
            phaseProperties.lookup<word>("systemType");
    }

    if (phaseCompressibleSystemType == word::null)
    {
        wordList phases
        (
            phaseProperties.lookupOrDefault("phases", wordList())
        );
        if (phases.size() == 2)
        {
            phaseCompressibleSystemType = "twoPhaseCompressibleSystem";
        }
        else if (phases.size() > 2)
        {
            phaseCompressibleSystemType = "multiphaseCompressibleSystem";
        }
        else
        {
            phaseCompressibleSystemType = "singlePhaseCompressibleSystem";
        }
    }
    Info<< "Selecting phaseCompressibleSystem: "
        << phaseCompressibleSystemType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(phaseCompressibleSystemType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseCompressibleSystem type "
            << phaseCompressibleSystemType << endl << endl
            << "Valid phaseCompressibleSystem types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}


// ************************************************************************* //
