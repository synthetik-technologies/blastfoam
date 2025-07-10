/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "polyRefinementConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "polyMeshPolyRefiner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(polyRefinementConstraint);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        polyRefinementConstraint,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyRefinementConstraint::polyRefinementConstraint
(
    const dictionary& constraintsDict,
    const word& modelType
)
:
    decompositionConstraint(constraintsDict, typeName)
{}


Foam::polyRefinementConstraint::polyRefinementConstraint()
:
    decompositionConstraint(dictionary(), typeName)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::polyRefinementConstraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    autoPtr<polyMeshPolyRefiner> storagePtr;
    polyMeshPolyRefiner* refPtr = nullptr;

    if (mesh.foundObject<polyMeshPolyRefiner>(polyMeshRefiner::typeName))
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : found polyMeshPolyRefiner" << endl;
        }
        refPtr = &mesh.lookupObjectRef<polyMeshPolyRefiner>
        (
            polyMeshRefiner::typeName
        );
    }
    else
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : reading polyMeshPolyRefiner from time "
                << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new polyMeshPolyRefiner(const_cast<polyMesh&>(mesh))
        );
    }

    polyMeshPolyRefiner& ref =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    // refinement itself implements decompositionConstraint
    ref.add
    (
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections
    );
}


void Foam::polyRefinementConstraint::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    autoPtr<polyMeshPolyRefiner> storagePtr;
    polyMeshPolyRefiner* refPtr = nullptr;

    if (mesh.foundObject<polyMeshPolyRefiner>(polyMeshRefiner::typeName))
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : found polyMeshPolyRefiner" << endl;
        }
        refPtr = &mesh.lookupObjectRef<polyMeshPolyRefiner>
        (
            polyMeshRefiner::typeName
        );
    }
    else
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : reading polyMeshPolyRefiner from time "
                << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new polyMeshPolyRefiner(const_cast<polyMesh&>(mesh))
        );
    }

    polyMeshPolyRefiner& ref =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    // refinement itself implements decompositionConstraint
    ref.apply
    (
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections,
        decomposition
    );
}


// ************************************************************************* //
