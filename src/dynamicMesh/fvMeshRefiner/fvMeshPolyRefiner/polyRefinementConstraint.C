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

#include "polyRefinementConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "fvMeshPolyRefiner.H"

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
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : setting constraints to preserve refinement history"
            << endl;
    }
}


Foam::polyRefinementConstraint::polyRefinementConstraint()
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

void Foam::polyRefinementConstraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    autoPtr<const fvMeshPolyRefiner> storagePtr;
    const fvMeshPolyRefiner* refPtr = nullptr;

    if (mesh.foundObject<fvMeshPolyRefiner>(fvMeshRefiner::typeName))
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : found fvMeshPolyRefiner" << endl;
        }
        refPtr = &mesh.lookupObject<fvMeshPolyRefiner>
        (
            fvMeshRefiner::typeName
        );
    }
    else
    {
//         if (decompositionConstraint::debug)
        {
            Info<< type() << " : reading fvMeshPolyRefiner from time "
                << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new fvMeshPolyRefiner
            (
                dynamicCast<fvMesh&>
                (
                    const_cast<polyMesh&>(mesh)
                )
            )
        );
    }

    const fvMeshPolyRefiner& ref =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    // refinement itself implements decompositionConstraint
    ref.refiner().add
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
    autoPtr<const fvMeshPolyRefiner> storagePtr;
    const fvMeshPolyRefiner* refPtr = nullptr;

    if (mesh.foundObject<fvMeshPolyRefiner>(fvMeshRefiner::typeName))
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : found fvMeshPolyRefiner" << endl;
        }
        refPtr = &mesh.lookupObject<fvMeshPolyRefiner>
        (
            fvMeshRefiner::typeName
        );
    }
    else
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : reading fvMeshPolyRefiner from time "
                << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new fvMeshPolyRefiner
            (
                dynamicCast<fvMesh&>
                (
                    const_cast<polyMesh&>(mesh)
                )
            )
        );
    }

    const fvMeshPolyRefiner& ref =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    // refinement itself implements decompositionConstraint
    ref.refiner().apply
    (
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections,
        decomposition
    );
}


// ************************************************************************* //
