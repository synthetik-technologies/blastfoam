/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-10-2022  Synthetik Applied Technologies : Collected functions
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

Namespace
    Foam::dynMeshTools

Description
    Collection of static functions to do various simple mesh related things.

SourceFiles
    dynMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef dynMeshTools_H
#define dynMeshTools_H

#include "meshTools.H"
#include "Pstream.H"
#include "transformField.H"
#include "Map.H"
#include "EdgeMap.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class IOobjectList;
class polyMesh;
class polyBoundaryMesh;
class coupledPolyPatch;
class polyTopoChange;
class faceZone;
class removeCells;

/*---------------------------------------------------------------------------*\
                        Namespace meshTools Declaration
\*---------------------------------------------------------------------------*/

namespace meshTools
{
    // Face information

        //- Set face information: patch, zone and zone flip for a face
        void getFaceInfo
        (
            const polyMesh& mesh,
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        );

        // Local topology modification functions (operate on cells/faces)

            //- Adds a face on top of existing faceI. Reverses if nessecary
            label addFace
            (
                polyTopoChange& meshMod,
                const polyMesh& mesh,
                const label faceI,
                const face& newFace,
                const label own,
                const label nei
            );

            //- Adds internal face from point. No checks on reversal
            label addInternalFace
            (
                polyTopoChange& meshMod,
                const polyMesh& mesh,
                const label meshFaceI,
                const label meshPointI,
                const face& newFace,
                const label own,
                const label nei
            );

            //- Modifies existing faceI for either new owner/neighbour or new
            //  face points. Reverses if nessecary
            void modifyFace
            (
                polyTopoChange& meshMod,
                const polyMesh& mesh,
                const label faceI,
                const face& newFace,
                const label own,
                const label nei
            );


        // Debug functions

            //- Check orientation of added internal face
            void checkInternalOrientation
            (
                const polyTopoChange& meshMod,
                const label cellI,
                const label faceI,
                const point& ownPt,
                const point& neiPt,
                const face& newFace
            );

            //- Check orientation of a new boundary face
            void checkBoundaryOrientation
            (
                const polyTopoChange& meshMod,
                const label cellI,
                const label faceI,
                const point& ownPt,
                const point& boundaryPt,
                const face& newFace
            );

            void checkFaceOrientation
            (
                const polyTopoChange& meshMod,
                const polyMesh& mesh,
                const label& faceI,
                const face& newFace
            );


        // Add an empty patch given a mesh, name, group, and dictionary
        label addPatch
        (
            fvMesh& mesh,
            const word& patchName,
            const wordList& groupNames,
            const dictionary& patchDict
        );

        label addPatch
        (
            fvMesh& mesh,
            const word& patchName,
            const word& group,
            const dictionary& patchDict
        );

        void modifyOrAddFace
        (
            polyTopoChange& meshMod,
            const face& f,
            const label facei,
            const label own,
            const bool flipFaceFlux,
            const label newPatchi,
            const label zoneID,
            const bool zoneFlip,

            PackedBoolList& modifiedFace
        );

        //- Create baffle faces for fZone faces. Usually newMasterPatches,
        //  newSlavePatches only size one but can be more for duplicate baffle
        //  sets
        label createBaffleFaces
        (
            const bool internalFacesOnly,
            const polyMesh& mesh,
            const faceZone& fZone,
            const labelList& newMasterPatches,
            const labelList& newSlavePatches,

            polyTopoChange& meshMod,
            PackedBoolList& modifiedFace
        );

        //- Create patch faces for fZone faces.
        label createPatchFaces
        (
            const bool internalFacesOnly,
            const polyMesh& mesh,
            const faceZone& fZone,
            const labelList& newPatches,

            polyTopoChange& meshMod,
            PackedBoolList& modifiedFace
        );

        //- Use the cells passed in to set the remove using remove cells
        label setRemoveCells
        (
            const polyMesh& mesh,
            const labelHashSet& selectedCells,
            const word& patchName,

            polyTopoChange& meshMod,
            const bool keepCells = false
        );

        //- Read and add fields to the database
        template<class FieldType>
        void readGeoFields(const fvMesh& mesh, const IOobjectList& objects);

        template<class FieldType>
        void readPointFields(const fvMesh& mesh, const IOobjectList& objects);

        //- Read and add all fields to the database
        void readAndStoreFields(const fvMesh& mesh);

        //- Remove empty patches
        void filterPatches
        (
            fvMesh& mesh,
            const HashSet<word>& addedPatchNames
        );
} // End namespace foamSyncTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dynMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
