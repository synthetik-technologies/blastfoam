/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
        void setFaceInfo
        (
            const polyMesh& mesh,
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
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
