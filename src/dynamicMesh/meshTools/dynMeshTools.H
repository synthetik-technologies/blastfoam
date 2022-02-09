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

/*
        //- Check whether uses couples.
        bool hasCouples(const polyBoundaryMesh&);

        //- Check for single transformation tensor only.
        void checkTransform(const coupledPolyPatch&, const bool);

        //- Apply separation to list. Either single vector or one vector
        //  per element.
        template <class T>
        void separateList(const vectorField&, UList<T>&);

        template <class T>
        void separateList(const vectorField&, Map<T>&);

        template <class T>
        void separateList(const vectorField&, EdgeMap<T>&);

        //- Combine value with existing value in map.
        template <class T, class CombineOp>
        void combine
        (
            Map<T>& pointValues,
            const CombineOp& cop,
            const label index,
            const T& val
        );

        //- Combine val with existing value at (implicit index) e.
        template <class T, class CombineOp>
        void combine
        (
            EdgeMap<T>& edgeValues,
            const CombineOp& cop,
            const edge& index,
            const T& val
        );

    // Static data members

        //- Synchronize values on all mesh points.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T, class CombineOp>
        void syncPointList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const T& nullValue,
            const bool applySeparation
        );

        //- Synchronize values on selected mesh points.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T, class CombineOp>
        void syncPointList
        (
            const polyMesh&,
            const labelList& meshPoints,
            UList<T>&,
            const CombineOp& bop,
            const T& nullValue,
            const bool applySeparation
        );

        //- Synchronize values on all mesh edges.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T, class CombineOp>
        void syncEdgeList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const T& nullValue,
            const bool applySeparation
        );

        //- Synchronize values on boundary faces only.
        //  Optionally applies rotation tensor for non-parallel cyclics
        //  (but not separation!)
        template <class T, class CombineOp>
        void syncBoundaryFaceList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const bool applySeparation
        );

        //- Synchronize values on all mesh faces.
        //  Optionally applies rotation tensor for non-parallel cyclics
        //  (but not separation!)
        template <class T, class CombineOp>
        void syncFaceList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const bool applySeparation
        );

        //- Swap coupled face values.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T>
        void swapBoundaryFaceList
        (
            const polyMesh&,
            UList<T>&,
            const bool applySeparation
        );

        //- Swap coupled face values.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T>
        void swapFaceList
        (
            const polyMesh&,
            UList<T>&,
            const bool applySeparation
        );


        // Sparse versions

            //- Synchronize values on selected points.
            //  Applies rotation and optionally separation for parallel
            //  cyclics.
            template <class T, class CombineOp>
            void syncPointMap
            (
                const polyMesh&,
                Map<T>& pointValues,
                const CombineOp& cop,
                const bool applySeparation
            );

            //- Synchronize values on selected edges. Edges are represented
            //  by the two vertices that make it up so global edges never get
            //  constructed.
            //  Applies rotation and optionally separation for parallel
            //  cyclics.
            template <class T, class CombineOp>
            void syncEdgeMap
            (
                const polyMesh&,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop,
                const bool applySeparation
            );


        // PackedList versions

            template <unsigned nBits, class CombineOp>
            void syncFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues,
                const CombineOp& cop
            );

            template <unsigned nBits>
            void swapFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues
            );

            template <unsigned nBits, class CombineOp>
            void syncPointList
            (
                const polyMesh& mesh,
                PackedList<nBits>& pointValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

            template <unsigned nBits, class CombineOp>
            void syncEdgeList
            (
                const polyMesh& mesh,
                PackedList<nBits>& edgeValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );


        // Other

            //- Get per point whether is it master (of a coupled set of points)
            PackedBoolList getMasterPoints(const polyMesh&);

            //- Get per edge whether is it master (of a coupled set of edges)
            PackedBoolList getMasterEdges(const polyMesh&);

            //- Get per face whether is it master (of a coupled set of faces)
            PackedBoolList getMasterFaces(const polyMesh&);

*/

} // End namespace foamSyncTools


/*
template<>
void foamSyncTools::separateList(const vectorField&, UList<vector>&);

template<>
void foamSyncTools::separateList(const vectorField&, Map<vector>&);

template<>
void foamSyncTools::separateList(const vectorField&, EdgeMap<vector>&);
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "dynMeshToolsTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
