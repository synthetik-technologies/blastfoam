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

#include "coupledCompressibleSystem.H"

#include "coupledSinglePhaseCompressibleSystem.H"
#include "coupledMultiphaseCompressibleSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::compressibleSystem>
Foam::coupledCompressibleSystem::New
(
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

    wordList phases
    (
        physicalPropertiesDict.lookupOrDefault("phases", wordList())
    );
    if (phases.size() < 2)
    {
        return autoPtr<compressibleSystem>
        (
            new coupledSinglePhaseCompressibleSystem
            (
                physicalPropertiesDict,
                mesh
            )
        );
    }
    else
    {
        return autoPtr<compressibleSystem>
        (
            new coupledMultiphaseCompressibleSystem
            (
                physicalPropertiesDict,
                mesh
            )
        );
    }
}


// ************************************************************************* //
