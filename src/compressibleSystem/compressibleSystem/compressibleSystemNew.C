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

#include "singlePhaseCompressibleSystem.H"
#include "twoPhaseCompressibleSystem.H"
#include "multiphaseCompressibleSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::compressibleSystem> Foam::compressibleSystem::New
(
    const fvMesh& mesh
)
{
    word compressibleSystemType(word::null);

    // Create temporary phase properties to lookup type
    // not store in the database to remove possible conflict
    Info<< "Reading phaseProperties dictionary\n" << endl;
    IOdictionary physicalPropertiesDict
    (
        IOobject
        (
            physicalProperties::typeName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    wordList phases
    (
        physicalPropertiesDict.lookupOrDefault("phases", wordList())
    );

    word ext = word::null;
    if (physicalPropertiesDict.found("sigma"))
    {
        if (physicalPropertiesDict.lookupOrDefault("useInterface", true))
        {
            ext = "Interface";
        }
    }
    else if (physicalPropertiesDict.lookupOrDefault("useInterface", false))
    {
        ext = "Interface";
    }

    if (phases.size() < 2)
    {
        return New
        (
            physicalPropertiesDict,
            mesh,
            singlePhaseCompressibleSystem::typeName,
            singlePhaseConstructorTablePtr_
        );
    }
    else if (phases.size() > 2)
    {
        return New
        (
            physicalPropertiesDict,
            mesh,
            multiphaseCompressibleSystem::typeName + ext,
            multiphaseConstructorTablePtr_
        );
    }
    return New
    (
        physicalPropertiesDict,
        mesh,
        twoPhaseCompressibleSystem::typeName + ext,
        twoPhaseConstructorTablePtr_
    );
}


Foam::autoPtr<Foam::compressibleSystem> Foam::compressibleSystem::NewCoupled
(
    const fvMesh& mesh
)
{
    word compressibleSystemType(word::null);

    // Create temporary phase properties to lookup type
    // not store in the database to remove possible conflict
    Info<< "Reading physicalProperties dictionary\n" << endl;
    IOdictionary physicalPropertiesDict
    (
        IOobject
        (
            physicalProperties::typeName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    wordList phases
    (
        physicalPropertiesDict.lookupOrDefault("phases", wordList())
    );

    // word ext = word::null;
    // if (physicalPropertiesDict.found("sigma"))
    // {
    //     if (physicalPropertiesDict.lookupOrDefault("useInterface", true))
    //     {
    //         ext = "Interface";
    //     }
    // }
    // else if (physicalPropertiesDict.lookupOrDefault("useInterface", false))
    // {
    //     ext = "Interface";
    // }


    return New
    (
        physicalPropertiesDict,
        mesh,
        (
            phases.size() < 2
          ? singlePhaseCompressibleSystem::typeName
          : multiphaseCompressibleSystem::typeName
        ) + "Coupled",
        coupledConstructorTablePtr_
    );
}


Foam::autoPtr<Foam::compressibleSystem> Foam::compressibleSystem::New
(
    const word& type,
    const fvMesh& mesh
)
{
    // Create temporary phase properties to lookup type
    // not store in the database to remove possible conflict
    Info<< "Reading physicalProperties dictionary\n" << endl;
    IOdictionary physicalPropertiesDict
    (
        IOobject
        (
            physicalProperties::typeName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Info<< "Selecting " << type << " compressibleSystem" << endl;
    {
        typename singlePhaseConstructorTable::iterator cstrIter =
            singlePhaseConstructorTablePtr_->find(type);
        if (cstrIter != singlePhaseConstructorTablePtr_->cend())
        {
            return cstrIter()(physicalPropertiesDict, mesh);
        }
    }
    {
        typename twoPhaseConstructorTable::iterator cstrIter =
            twoPhaseConstructorTablePtr_->find(type);
        if (cstrIter != twoPhaseConstructorTablePtr_->cend())
        {
            return cstrIter()(physicalPropertiesDict, mesh);
        }
    }
    {
        typename multiphaseConstructorTable::iterator cstrIter =
            multiphaseConstructorTablePtr_->find(type);
        if (cstrIter != multiphaseConstructorTablePtr_->cend())
        {
            return cstrIter()(physicalPropertiesDict, mesh);
        }
    }


    FatalErrorInFunction
        << "Unknown compressibleSystem type " << type << endl << endl
        << "Valid compressibleSystem types are : " << endl
        << singlePhaseConstructorTablePtr_->sortedToc() << nl
        << twoPhaseConstructorTablePtr_->sortedToc() << nl
        << multiphaseConstructorTablePtr_->sortedToc() << nl
        << exit(FatalError);

    return singlePhaseConstructorTablePtr_->find(type)()
    (
        physicalPropertiesDict,
        mesh
    );
}

// ************************************************************************* //
