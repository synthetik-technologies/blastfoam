/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "hexRefRefinementHistoryConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "hexRefRefinementHistory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(hexRefRefinementHistoryConstraint);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        hexRefRefinementHistoryConstraint,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexRefRefinementHistoryConstraint::hexRefRefinementHistoryConstraint
(
    const dictionary& constraintsDict,
    const word& modelType
)
:
    decompositionConstraint(constraintsDict, typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : setting constraints to preserve refinement history"
            << endl;
    }
}


Foam::hexRefRefinementHistoryConstraint::hexRefRefinementHistoryConstraint()
:
    decompositionConstraint(dictionary(), typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : setting constraints to refinement history"
            << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hexRefRefinementHistoryConstraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    autoPtr<const hexRefRefinementHistory> storagePtr;
    hexRefRefinementHistory const* refPtr = nullptr;

    if (mesh.foundObject<hexRefRefinementHistory>("refinementHistory"))
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : found refinementHistory" << endl;
        }
        refPtr = &mesh.lookupObject<hexRefRefinementHistory>("refinementHistory");
    }
    else
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : reading refinementHistory from time "
                << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new hexRefRefinementHistory
            (
                IOobject
                (
                    "refinementHistory",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh.nCells()
            )
        );
    }

    const hexRefRefinementHistory& history =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    if (history.active())
    {
        // hexRefRefinementHistory itself implements decompositionConstraint
        history.add
        (
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections
        );
    }
}


void Foam::hexRefRefinementHistoryConstraint::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    autoPtr<const hexRefRefinementHistory> storagePtr;
    hexRefRefinementHistory const* refPtr = nullptr;

    if (mesh.foundObject<hexRefRefinementHistory>("refinementHistory"))
    {
        // if (decompositionConstraint::debug)
        //{
        //    Info<< type() << " : found refinementHistory" << endl;
        //}
        refPtr = &mesh.lookupObject<hexRefRefinementHistory>("refinementHistory");
    }
    else
    {
        // if (decompositionConstraint::debug)
        //{
        //    Info<< type() << " : reading refinementHistory from time "
        //        << mesh.facesInstance() << endl;
        //}
        storagePtr.reset
        (
            new hexRefRefinementHistory
            (
                IOobject
                (
                    "refinementHistory",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh.nCells()
            )
        );
    }

    const hexRefRefinementHistory& history =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    if (history.active())
    {
        // hexRefRefinementHistory itself implements decompositionConstraint
        history.apply
        (
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections,
            decomposition
        );
    }
}


// ************************************************************************* //
