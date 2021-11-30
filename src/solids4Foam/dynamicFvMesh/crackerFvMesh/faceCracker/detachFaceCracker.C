/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "faceCracker.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "processorPolyPatch.H"
#include "labelPair.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceCracker::detachFaceCracker
(
    polyTopoChange& ref
) const
{
    FatalErrorInFunction
        << "No longer implemented: use instead detachInternalFaces"
        << "and detachCoupledFaces" << abort(FatalError);
}


void Foam::faceCracker::detachInternalFaces
(
    polyTopoChange& ref
) const
{
    // Method
    // 1) Check if edges of the face to break are internal:
    //        Three edge types: standard edges, proc edges, global edges.
    // 2) Find added points
    //        Points that are only on a boundary edge
    // 3) Check if surrounding faces need to be modified
    //        First add all point-faces of points-to-add
    //        Then, starting at the master face (face-to-break), perform a
    //        cell-face-cell walk, through faces with a point-to-add, removing
    //        faces from the faces-to-modify list
    // 4) If faces were modified
    //        Add the points-to-add
    // 5) Modify surrounding faces
    // 6) Modify the master face to the crack patch, and add a slave face


    if (debug)
    {
        Pout<< nl << "detachInternalFaces" << nl << endl;
    }


    // 1) Check if edges of the face to break are internal:
    //        Three edge types: standard edges, proc edges, global edges.

    const polyMesh& mesh = topoChanger().mesh();
    const meshFaceZones& zoneMesh = mesh.faceZones();

    if (zoneMesh[crackZoneID_.index()].size() > 1)
    {
        FatalErrorInFunction
            << "Only one internal face can be broken at a time"
            << abort(FatalError);
    }

    {
        // In parallel, the face-to-break will be on just one processor;
        // But it may have edges at a processor boundary
        label faceToBreakID = -1;
        face faceToBreak(0);
        bool faceToBreakFlip = false;
        labelList curFaceEdges(0);
        const labelList& faceOwn = mesh.faceOwner();
        label faceCellID = -1;
        if (zoneMesh[crackZoneID_.index()].size())
        {
            faceToBreakID = zoneMesh[crackZoneID_.index()][0];
            faceToBreak = mesh.faces()[faceToBreakID];
            faceToBreakFlip = zoneMesh[crackZoneID_.index()].flipMap()[0];
            curFaceEdges = mesh.faceEdges()[faceToBreakID];
            faceCellID = faceOwn[faceToBreakID];

            //if (debug)
            {
                Pout<< "Breaking internal face : "
                    << mesh.faceCentres()[faceToBreakID] << endl;
            }
        }
        boolList edgeIsInternal(faceToBreak.nEdges(), true);
        const labelListList& edgeFaces = mesh.edgeFaces();
        const edgeList& edges = mesh.edges();
        const labelList& faceNei = mesh.faceNeighbour();
        const faceList& faces = mesh.faces();
        const labelListList& pointFaces = mesh.pointFaces();

        // Edges shared by two processors
        labelList procEdge(edgeIsInternal.size(), -1);
        labelList procEdgePatch(edgeIsInternal.size(), -1);
        label nProcEdges = 0;

        // Edges of processor boundaries shared by more than two processors (so
        // called global edges)
        const labelList& glEdges = mesh.globalData().sharedEdgeLabels();
        const labelList& glEdgeAddr = mesh.globalData().sharedEdgeAddr();
        const labelList& glPoints = mesh.globalData().sharedPointLabels();
        const labelList& glPointAddr = mesh.globalData().sharedPointAddr();
        labelHashSet sharedEdgeSet;

        // Check internal edges fully on the current processor
        forAll(curFaceEdges, eI)
        {
            const label edgeID = curFaceEdges[eI];
            const labelList& curEdgeFaces = edgeFaces[edgeID];

            forAll (curEdgeFaces, fI)
            {
                const label faceID = curEdgeFaces[fI];

                if (!mesh.isInternalFace(faceID))
                {
                    const label patchID =
                        mesh.boundaryMesh().whichPatch(faceID);

                    const polyPatch& ppatch = mesh.boundaryMesh()[patchID];

                    if (!isA<processorPolyPatch>(ppatch))
                    {
                        // The edge belongs to a boundary face
                        edgeIsInternal[eI] = false;
                    }
                    else
                    {
                        // Check if the face is a global edge i.e. if it is
                        // shared by more than two processors
                        bool glEdge = (findIndex(glEdges, edgeID) != -1);

                        if (glEdge)
                        {
                            // This is a global edge and requires special
                            // treatment
                            if (!sharedEdgeSet.found(edgeID))
                            {
                                sharedEdgeSet.insert(edgeID);
                            }
                        }
                        else
                        {
                            // This edge is on a processor boundary
                            procEdge[eI] =
                                findIndex(ppatch.meshEdges(), edgeID);
                            procEdgePatch[eI] = patchID;

                            nProcEdges++;
                        }
                    }
                }
            }
        }

        reduce(nProcEdges, maxOp<label>());

        // Check processor edges
        // An edge is only internal if it is internal on both processors

        if (nProcEdges)
        {
            labelList receivedProcEdge(0);
            labelList receivedProcEdgePatch(0);
            boolList receivedProcEdgeIsInternal(0);

            // Send edges to the neighbour processor to be checked
            forAll(mesh.boundaryMesh(), patchI)
            {
                const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

                if (isA<processorPolyPatch>(ppatch))
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(ppatch);

                    // Find edges just on this patch

                    labelPairHashSet curPatchEdgesSet;

                    forAll(procEdge, eI)
                    {
                        if (procEdgePatch[eI] == patchI)
                        {
                            label edgeIsInt = 0;
                            if (edgeIsInternal[eI])
                            {
                                edgeIsInt = 1;
                            }

                            labelPair curPair(procEdge[eI], edgeIsInt);

                            curPatchEdgesSet.insert(curPair, nil());
                        }
                    }

                    List<labelPair> curPatchProcEdge = curPatchEdgesSet.toc();

                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo()
                    );

                    toProc
                        << curPatchProcEdge;
                }
            }


            // Receive edges from the neighbour processor
            forAll(mesh.boundaryMesh(), patchI)
            {
                const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

                if (isA<processorPolyPatch>(ppatch))
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(ppatch);

                    // Find edges just on this patch

                    List<labelPair> receivedPatchProcEdge;

                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo()
                    );

                    fromProc
                        >> receivedPatchProcEdge;

                    // Convert receovedProcEdges into the current patch
                    // addressing
                    forAll(receivedPatchProcEdge, eI)
                    {
                        labelPair& curPair = receivedPatchProcEdge[eI];
                        const label edgeID = curPair.first();

                        const label curSize = receivedProcEdge.size();
                        receivedProcEdge.setSize(curSize + 1);
                        receivedProcEdgeIsInternal.setSize(curSize + 1);
                        receivedProcEdgePatch.setSize(curSize + 1);

                        receivedProcEdge[curSize] =
                            findIndex
                            (
                                procPatch.nbrEdges(),
                                edgeID
                            );

                        receivedProcEdgeIsInternal[curSize] =
                            bool(curPair.second() == 1);

                        receivedProcEdgePatch[curSize] = patchI;
                    }
                }
            }

            // Now we check if the edges received from the neighbour processor
            // are internal

            forAll(receivedProcEdge, eI)
            {
                const label edgeID = receivedProcEdge[eI];
                const label edgePatchID = receivedProcEdgePatch[eI];

                const labelList& curMeshEdges =
                    mesh.boundaryMesh()[edgePatchID].meshEdges();

                const labelList& curFaces = edgeFaces[curMeshEdges[edgeID]];

                forAll(curFaces, faceI)
                {
                    if (!mesh.isInternalFace(curFaces[faceI]))
                    {
                        const label patchID =
                            mesh.boundaryMesh().whichPatch(curFaces[faceI]);

                        if
                        (
                           !isA<processorPolyPatch>
                            (
                                mesh.boundaryMesh()[patchID]
                            )
                        )
                        {
                            // The edge belongs to a boundary face
                            receivedProcEdgeIsInternal[eI] = false;
                            break;
                        }
                    }
                }
            }

            // Now send the receivedProcEdgeIsInternal list back to the
            // neighbour processor

            forAll(mesh.boundaryMesh(), patchI)
            {
                const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

                if (isA<processorPolyPatch>(ppatch))
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(ppatch);

                    labelPairHashSet curPatchEdgesSet;

                    forAll(receivedProcEdge, eI)
                    {
                        if (receivedProcEdgePatch[eI] == patchI)
                        {
                            label edgeIsInt = 0;
                            if (receivedProcEdgeIsInternal[eI])
                            {
                                edgeIsInt = 1;
                            }

                            labelPair curPair(receivedProcEdge[eI], edgeIsInt);

                            curPatchEdgesSet.insert(curPair, nil());
                        }
                    }

                    List<labelPair> curPatchProcEdge = curPatchEdgesSet.toc();

                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo()
                    );

                    toProc
                        << curPatchProcEdge << endl;
                }
            }


            // Receive edges and check if internal
            // Reset received fields
            receivedProcEdge.setSize(0);
            receivedProcEdgePatch.setSize(0);
            receivedProcEdgeIsInternal.setSize(0);

            forAll(mesh.boundaryMesh(), patchI)
            {
                const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

                if (isA<processorPolyPatch>(ppatch))
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(ppatch);

                    List<labelPair> receivedPatchProcEdge;

                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo()
                    );

                    fromProc
                        >> receivedPatchProcEdge;

                    // Convert receovedProcEdges into the current patch
                    // addressing
                    forAll(receivedPatchProcEdge, eI)
                    {
                        labelPair& curPair = receivedPatchProcEdge[eI];
                        const label edgeID = curPair.first();

                        const label curSize = receivedProcEdge.size();
                        receivedProcEdge.setSize(curSize + 1);
                        receivedProcEdgeIsInternal.setSize(curSize + 1);
                        receivedProcEdgePatch.setSize(curSize + 1);

                        receivedProcEdge[curSize] =
                            findIndex
                            (
                                procPatch.nbrEdges(),
                                edgeID
                            );

                        receivedProcEdgeIsInternal[curSize] =
                            bool(curPair.second() == 1);

                        receivedProcEdgePatch[curSize] = patchI;
                    }
                }
            }

            // Sync the edgeIsInternal list: an edge is internal if all
            // processors say it is internal
            forAll(procEdge, eI)
            {
                const label edgeID = procEdge[eI];

                // Check if this edge is in the received list
                forAll(receivedProcEdge, reI)
                {
                    const label rEdgeID = receivedProcEdge[reI];

                    if (edgeID == rEdgeID)
                    {
                        if (!receivedProcEdgeIsInternal[reI])
                        {
                            edgeIsInternal[eI] = false;
                        }
                    }
                }
            }
        }


        // Check global edges
        {
            const labelList sharedEdges = sharedEdgeSet.toc();

            // Note: we use a labelList instead of a boolList because there is
            // no reduce(..., orOp<boolList>(...)) function
            labelList checkSharedEdge(mesh.globalData().nGlobalEdges(), 0);

            forAll(sharedEdges, eI)
            {
                const label edgeID = sharedEdges[eI];

                const label glEdgeID = findIndex(glEdges, edgeID);

                if (glEdgeID != -1)
                {
                    checkSharedEdge[glEdgeAddr[glEdgeID]] = 1;
                }
            }

            reduce(checkSharedEdge, maxOp<labelList>());

            labelList internalSharedEdge(mesh.globalData().nGlobalEdges(), 1);

            forAll(checkSharedEdge, glEdgeI)
            {
                if (checkSharedEdge[glEdgeI])
                {
                    const label sharedEdgeIndex =
                        findIndex(glEdgeAddr, glEdgeI);

                    if (sharedEdgeIndex != -1)
                    {
                        const label curSharedEdge = glEdges[sharedEdgeIndex];

                        const labelList& curFaces = edgeFaces[curSharedEdge];

                        // Check if the edge is internal on this processor
                        forAll(curFaces, fI)
                        {
                            const label faceID = curFaces[fI];

                            if (!mesh.isInternalFace(faceID))
                            {
                                const label patchID =
                                    mesh.boundaryMesh().whichPatch(faceID);

                                if
                                (
                                   !isA<processorPolyPatch>
                                    (
                                        mesh.boundaryMesh()[patchID]
                                    )
                                )
                                {
                                    // The edge belongs to a boundary face
                                    internalSharedEdge[glEdgeI] = 0;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            reduce(internalSharedEdge, minOp<labelList>());

            forAll(sharedEdges, eI)
            {
                const label edgeID = sharedEdges[eI];

                const label glEdgeID = findIndex(glEdges, edgeID);

                if (glEdgeID != -1)
                {
                    // Update local edgeIsInternal list
                    forAll(curFaceEdges, faceToBreakEdgeI)
                    {
                        if (curFaceEdges[faceToBreakEdgeI] == edgeID)
                        {
                            if (edgeIsInternal[faceToBreakEdgeI])
                            {
                                edgeIsInternal[faceToBreakEdgeI] =
                                    internalSharedEdge[glEdgeAddr[glEdgeID]];
                                break;
                            }
                        }
                    }
                }
            }
        }

        // Print out which edges are internal
        if (debug)
        {
            const pointField& points = mesh.points();

            Pout<< nl << "internal edges:" << endl;

            forAll(edgeIsInternal, eI)
            {
                Pout<< "    " << edges[curFaceEdges[eI]].centre(points)
                    << " " << edgeIsInternal[eI] << endl;
            }
        }


        // 2) Find added points
        //        Points that are only on a boundary edge

        labelHashSet pointsToAddSet(faceToBreak.size());

        // Initially add all points of the face-to-break, then we will remove
        // those on an internal edge
        forAll(faceToBreak, pI)
        {
            pointsToAddSet.insert(faceToBreak[pI]);
        }

        forAll(edgeIsInternal, eI)
        {
            if (edgeIsInternal[eI])
            {
                const label edgeID = curFaceEdges[eI];
                const edge& curEdge = edges[edgeID];

                // Remove points from list if found

                if (pointsToAddSet.found(curEdge.start()))
                {
                    pointsToAddSet.erase(curEdge.start());
                }

                if (pointsToAddSet.found(curEdge.end()))
                {
                    pointsToAddSet.erase(curEdge.end());
                }
            }
        }


        // Send split processor boundary points to the neighbour processor so it
        // knows to split them too

        const labelList pointsToAddNoSync = pointsToAddSet.toc();

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(ppatch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(ppatch);

                labelHashSet curPatchPointsToAddSet;

                forAll(pointsToAddNoSync, pI)
                {
                    const label pointID = pointsToAddNoSync[pI];

                    const label localPointID = ppatch.whichPoint(pointID);

                    if (localPointID != -1)
                    {
                        curPatchPointsToAddSet.insert(localPointID);
                    }
                }

                labelList curPatchPointsToAdd = curPatchPointsToAddSet.toc();

                OPstream toProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                );

                toProc
                    << curPatchPointsToAdd << endl;
            }
        }


        // Receive pointsToAdd from the neighbour processor
        labelList pointsToAddFromNeiProc(0);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(ppatch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(ppatch);

                labelList receivedPointsToAddFromNeiProc;

                IPstream fromProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                );

                fromProc
                    >> receivedPointsToAddFromNeiProc;

                // Convert receivedPointsToAddFromNeiProc into the current patch
                // addressing and convert to mesh indexing

                const labelList& pMeshPoints = ppatch.meshPoints();

                forAll(receivedPointsToAddFromNeiProc, pI)
                {
                    const label pointID = receivedPointsToAddFromNeiProc[pI];

                    const label curSize = pointsToAddFromNeiProc.size();
                    pointsToAddFromNeiProc.setSize(curSize + 1);

                    const label curProcPointID =
                        pMeshPoints
                        [
                            findIndex
                            (
                                procPatch.nbrPoints(),
                                pointID
                            )
                        ];

                    pointsToAddFromNeiProc[curSize] = curProcPointID;

                    // Insert into the pointsToAdd list
                    if (!pointsToAddSet.found(curProcPointID))
                    {
                        pointsToAddSet.insert(curProcPointID);
                    }
                }
            }
        }


        // Add global points to the pointsToAddSet
        labelHashSet sharedPointsToAddSet;
        {
            const labelList pointsToAddBeforeGlobalCheck = pointsToAddSet.toc();

            forAll(pointsToAddBeforeGlobalCheck, pI)
            {
                const label pointID = pointsToAddBeforeGlobalCheck[pI];
                const label glPointID = findIndex(glPoints, pointID);

                if (glPointID != -1)
                {
                    sharedPointsToAddSet.insert(pointID);
                }
            }
        }

        const labelList sharedPointsToAdd = sharedPointsToAddSet.toc();

        labelList checkSharedPoint(mesh.globalData().nGlobalPoints(), 0);

        forAll(sharedPointsToAdd, pI)
        {
            const label pointID = sharedPointsToAdd[pI];

            const label glPointID = findIndex(glPoints, pointID);

            if (glPointID != -1)
            {
                checkSharedPoint[glPointAddr[glPointID]] = 1;
            }
        }

        reduce(checkSharedPoint, maxOp<labelList>());

        forAll(checkSharedPoint, glPointI)
        {
            if (checkSharedPoint[glPointI])
            {
                const label sharedPointIndex = findIndex(glPointAddr, glPointI);

                if (sharedPointIndex != -1)
                {
                    const label pointID = glPoints[sharedPointIndex];

                    // Add point to be split
                    if (!pointsToAddSet.found(pointID))
                    {
                        pointsToAddSet.insert(pointID);
                    }
                }
            }
        }


        //const labelList pointsToAdd = pointsToAddSet.toc();
        labelList pointsToAdd = pointsToAddSet.toc();

        // Print points-to-add
        if (debug)
        {
            const pointField& points = mesh.points();

            Pout<< nl << "pointsToAdd: " << endl;

            forAll(pointsToAdd, pI)
            {
                Pout<< "    " << points[pointsToAdd[pI]] << endl;
            }
        }



       // 3) Check if surrounding faces need to be modified
        //        First add all point-faces of points-to-add
        //        Then, starting at the master face (face-to-break), perform a
        //        cell-face-cell walk, through faces with a point-to-add,
        //        removing faces from the faces-to-modify list

        // A processor may have a point to add, but if it does not need to
        // modify faces then we do not need to add a point
        // Try find a boundary face attached to the points-to-add; if we cannot
        // find a "master face" then we will not add points on this proc
        if (faceCellID == -1 && pointsToAdd.size())
        {
            // This proc does not have a face to break but it does need to
            // split an edge
            // Pick a boundary face to be the master face
            label masterFaceID = -1;
            forAll(pointsToAdd, pI)
            {
                const label pointID = pointsToAdd[pI];
                const labelList& curPointFaces = pointFaces[pointID];

                forAll(curPointFaces, fI)
                {
                    const label faceID = curPointFaces[fI];

                    if (!mesh.isInternalFace(faceID))
                    {
                        const label patchID =
                            mesh.boundaryMesh().whichPatch(faceID);

                        if
                        (
                            !isA<processorPolyPatch>
                            (
                                mesh.boundaryMesh()[patchID]
                            )
                        )
                        {
                            // This will be the master face
                            masterFaceID = faceID;
                            break;
                        }
                    }
                }

                if (masterFaceID != -1)
                {
                    break;
                }
            }

            if (masterFaceID != -1)
            {
                faceCellID = faceOwn[masterFaceID];
            }
        }


        labelHashSet facesToModifySet(20);

        if (faceCellID != -1)
        {
            forAll(pointsToAdd, pI)
            {
                const labelList& curPointFaces = pointFaces[pointsToAdd[pI]];

                forAll(curPointFaces, pfI)
                {
                    const label pointFaceID = curPointFaces[pfI];

                    if (!facesToModifySet.found(pointFaceID))
                    {
                        facesToModifySet.insert(pointFaceID);
                    }
                }
            }

            // Starting at the master faceCell (owner of the face to break),
            // check all faces of the cell; if the face is in the facesToModify
            // list then remove it; add the neighbour cell to be checked next.
            // We will continue this cell-face-cell walk until there are no more
            // cells to check.

            SLList<label> cellsToCheck;
            cellsToCheck.append(faceCellID);
            labelHashSet checkedFaces(20);

            // Add faceToBreak to checked faces to stop cell-walk through this
            // face and remove it from the facesToModify
            checkedFaces.insert(faceToBreakID);
            facesToModifySet.erase(faceToBreakID);

            const cellList& cells = mesh.cells();

            // Cell-face-cell walk
            do
            {
                const label cellID = cellsToCheck.first();
                const labelList& curCellFaces = cells[cellID];

                forAll(curCellFaces, fI)
                {
                    const label faceID = curCellFaces[fI];

                    if (!checkedFaces.found(faceID))
                    {
                        checkedFaces.insert(faceID);

                        // Remove from the facesToModifySet if found
                        if (facesToModifySet.found(faceID))
                        {
                            facesToModifySet.erase(faceID);

                            // Add neighbour cell to cellsToCheck
                            if (mesh.isInternalFace(faceID))
                            {
                                label neiCellID = faceNei[faceID];
                                if (neiCellID == cellID)
                                {
                                    neiCellID = faceOwn[faceID];
                                }

                                cellsToCheck.append(neiCellID);
                            }
                        }
                    }
                }

                // Remove current cell from cellsToCheck
                cellsToCheck.removeHead();
            }
            while (cellsToCheck.size());
        }

        labelList addedPoints(pointsToAdd.size(), -1);

        // Print faces to be modified
        if (facesToModifySet.size() == 0)
        {
            addedPoints.setSize(0);
            pointsToAdd.setSize(0);
        }
        else //if (facesToModifySet.size())
        {
            const labelList facesToModify = facesToModifySet.toc();

            if (debug)
            {
                Pout<< nl << "facesToModify: " << endl;

                forAll(facesToModify, fI)
                {
                    Pout<< "Modify face "
                        << mesh.faceCentres()[facesToModify[fI]]
                        << endl;
                }
            }


            // 4) If faces were modified
            //        Add the points-to-add

            const pointField& points = mesh.points();

            forAll(pointsToAdd, pI)
            {
                const label pointID = pointsToAdd[pI];

                addedPoints[pI] =
                    ref.setAction
                    (
                        polyAddPoint
                        (
                            points[pointID],           // point
                            pointID,                   // master point
                            -1,                        // zone ID
                            true                       // supports a cell
                        )
                    );

                if (debug)
                {
                    Pout<< "Adding point " << points[pointID] << endl;
                }
            }


            // 5) Modify surrounding faces

            forAll(facesToModify, fI)
            {
                const label faceID = facesToModify[fI];

                const face& oldFace = faces[faceID];
                face newFace = oldFace;

                forAll (oldFace, pI)
                {
                    const label oldPointID = oldFace[pI];

                    // Check if the point has been split
                    forAll(pointsToAdd, ptoaI)
                    {
                        const label pointToAddID = pointsToAdd[ptoaI];

                        if (pointToAddID == oldPointID)
                        {
                            newFace[pI] = addedPoints[ptoaI];
                            break;
                        }
                    }
                }

                if (mesh.isInternalFace(faceID))
                {
                    ref.setAction
                    (
                        polyModifyFace
                        (
                            newFace,                    // face
                            faceID,                     // master face
                            faceOwn[faceID],            // owner
                            faceNei[faceID],            // neighbour
                            false,                      // flip flux
                            -1,                         // patch for face
                            false,                      // remove from zone
                            -1,                         // zone for face
                            false                       // face zone flip
                        )
                    );
                }
                else
                {
                    ref.setAction
                    (
                        polyModifyFace
                        (
                            newFace,                     // face
                            faceID,                      // master face
                            faceOwn[faceID],             // owner
                            -1,                          // neighbour
                            false,                       // flip flux
                            mesh.boundaryMesh().whichPatch(faceID), // patch
                            false,                        // remove from zone
                            -1,                           // zone for face
                            false                         // face zone flip
                        )
                    );
                }
            }
        }


        // 6) Move the master face to the crack patch and add the slave face

        // Move master face to the crack patch
        // This is only performed by the processor with the internal face to
        // break

        if (faceToBreakID != -1)
        {
            const label faceCellID = faceOwn[faceToBreakID];
            const label faceNeiCellID = faceNei[faceToBreakID];

            if (faceToBreakFlip)
            {
                // Face needs to be flipped for the master patch
                ref.setAction
                (
                    polyModifyFace
                    (
                        faceToBreak.reverseFace(), // modified face
                        faceToBreakID,           // label of face being modified
                        faceNeiCellID,            // owner
                        -1,                       // neighbour
                        true,                     // face flip
                        crackPatchID_.index(),    // patch for face
                        false,                    // remove from zone
                        crackZoneID_.index(),     // zone for face
                        !faceToBreakFlip          // face flip in zone
                    )
                );
            }
            else
            {
                // No flip
                ref.setAction
                (
                    polyModifyFace
                    (
                        faceToBreak,              // modified face
                        faceToBreakID,           // label of face being modified
                        faceCellID,               // owner
                        -1,                       // neighbour
                        false,                    // face flip
                        crackPatchID_.index(),    // patch for face
                        false,                    // remove from zone
                        crackZoneID_.index(),     // zone for face
                        faceToBreakFlip           // face flip in zone
                    )
                );
            }

            // Add the slave face

            // Build the face for the slave patch by renumbering
            const face oldFace = faceToBreak.reverseFace();
            face newFace = oldFace;

            forAll (oldFace, pI)
            {
                const label oldPointID = oldFace[pI];

                // Check if the point has been split
                forAll(pointsToAdd, ptoaI)
                {
                    const label pointToAddID = pointsToAdd[ptoaI];

                    if (pointToAddID == oldPointID)
                    {
                        newFace[pI] = addedPoints[ptoaI];
                        break;
                    }
                }
            }

            if (min(newFace) == -1)
            {
                FatalErrorIn("detachInternalFace()")
                    << "newFace: " << newFace << nl
                    << "oldFace.reverseFace(): " << oldFace.reverseFace() << nl
                    << "addedPoints: " << addedPoints << nl
                    << "pointsToAdd: " << pointsToAdd << nl
                    << abort(FatalError);
            }

            if (faceToBreakFlip)
            {
                // Add slave face
                ref.setAction
                (
                    polyAddFace
                    (
                        newFace,                        // face
                        faceCellID,                     // owner
                        -1,                             // neighbour
                        -1,                             // master point
                        -1,                             // master edge
                        faceToBreakID,                  // master face
                        false,                          // flip flux
                        crackPatchID_.index(),          // new patch index
                        -1,                             // zone for face
                        false                           // zone flip
                    )
                );
            }
            else
            {
                // Add renumbered face into the slave patch
                ref.setAction
                (
                    polyAddFace
                    (
                        newFace,                        // face
                        faceNeiCellID,                  // owner
                        -1,                             // neighbour
                        -1,                             // master point
                        -1,                             // master edge
                        faceToBreakID,                  // master face
                        true,                           // flip flux
                        crackPatchID_.index(),          // new patch index
                        -1,                             // zone for face
                        false                           // face flip in zone
                    )
                );
            }
        }
    }
}


void Foam::faceCracker::detachCoupledFaces
(
    polyTopoChange& ref
) const
{
    // Method
    // 1) Check if edges of the face to break are internal:
    //        Three edge types: standard edges, proc edges, global edges.
    // 2) Find added points
    //        Points that are only on a boundary edge
    // 3) Check if surrounding faces need to be modified
    //        First add all point-faces of points-to-add
    //        Then, starting at the master face (face-to-break), perform a
    //        cell-face-cell walk, through faces with a point-to-add, removing
    //        faces from the faces-to-modify list
    // 4) If faces were modified
    //        Add the points-to-add
    // 5) Modify surrounding faces
    // 6) Move the coupled face to the crack patch

    if (debug)
    {
        Pout<< nl << "detachCoupledFaces" << nl << endl;
    }

    // 1) Check if edges of the face to break are internal:
    //        Three edge types: standard edges, proc edges, global edges.

    if (coupledFacesToBreak_.size() > 1)
    {
        FatalErrorIn("faceCracker::detachCoupledFaces()")
            << "Only one coupled face can be broken at a time"
            << abort(FatalError);
    }

    {
        const polyMesh& mesh = topoChanger().mesh();

        label faceToBreakID = -1;
        face faceToBreak(0);
        labelList curFaceEdges(0);
        label procPatchID = -1;
        label faceCellID = -1;
        const labelList& faceOwn = mesh.faceOwner();
        if (coupledFacesToBreak_.size())
        {
            faceToBreakID = coupledFacesToBreak_[0];
            faceToBreak = mesh.faces()[faceToBreakID];
            curFaceEdges = mesh.faceEdges()[faceToBreakID];
            procPatchID = mesh.boundaryMesh().whichPatch(faceToBreakID);
            faceCellID = faceOwn[faceToBreakID];

            //if (debug)
            {
                Pout<< "Breaking coupled face : "
                    << mesh.faceCentres()[faceToBreakID] << endl;
            }
        }
        boolList edgeIsInternal(faceToBreak.nEdges(), true);
        const labelListList& edgeFaces = mesh.edgeFaces();
        const labelListList& pointFaces = mesh.pointFaces();
        const faceList& faces = mesh.faces();
        const edgeList& edges = mesh.edges();
        const labelList& faceNei = mesh.faceNeighbour();

        // Edges shared by two processors
        labelList procEdge(edgeIsInternal.size(), -1);

        // Edges of processor boundaries shared by more than two processors
        const labelList& glEdges = mesh.globalData().sharedEdgeLabels();
        const labelList& glEdgeAddr = mesh.globalData().sharedEdgeAddr();
        const labelList& glPoints = mesh.globalData().sharedPointLabels();
        const labelList& glPointAddr = mesh.globalData().sharedPointAddr();
        labelHashSet sharedEdgeSet;

        // Check if edges are internal on the current processor
        forAll(curFaceEdges, eI)
        {
            const label edgeID = curFaceEdges[eI];
            const labelList& curEdgeFaces = edgeFaces[edgeID];

            forAll(curEdgeFaces, fI)
            {
                const label faceID = curEdgeFaces[fI];

                if (!mesh.isInternalFace(faceID))
                {
                    const label patchID =
                        mesh.boundaryMesh().whichPatch(faceID);

                    const polyPatch& ppatch = mesh.boundaryMesh()[patchID];

                    if (!isA<processorPolyPatch>(ppatch))
                    {
                        // The edge belongs to a boundary face
                        edgeIsInternal[eI] = false;
                    }
                    else
                    {
                        // Check if the face is a global edge i.e. if it is
                        // shared by more than two processors
                        bool glEdge = (findIndex(glEdges, edgeID) != -1);

                        if (glEdge)
                        {
                            // This is a global edge and requires special
                            // treatment
                            if (!sharedEdgeSet.found(edgeID))
                            {
                                sharedEdgeSet.insert(edgeID);
                            }
                        }
                        else
                        {
                            // This edge is on a processor boundary
                            procEdge[eI] =
                                findIndex(ppatch.meshEdges(), edgeID);
                        }
                    }
                }
            }
        }

        // Check processor edges
        // An edge is only internal if it is internal on both processors
        if (faceToBreakID != -1)
        {
            const polyPatch& ppatch = mesh.boundaryMesh()[procPatchID];
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(ppatch);
            labelList receivedProcEdge(procEdge.size(), -1);
            boolList receivedProcEdgeIsInternal(edgeIsInternal.size(), true);

            // Send edges to the neighbour processor
            {
                OPstream toProc(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                toProc
                    << procEdge << edgeIsInternal;
            }

            // Receive edges from the neighbour processor
            {
                IPstream fromProc(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                fromProc
                    >> receivedProcEdge >> receivedProcEdgeIsInternal;

                // Convert receivedProcEdges into the current patch addressing
                // The points/edges on the neighbour processor patch are in the
                // same order as on the current processor patch
                forAll(receivedProcEdge, eI)
                {
                    if (receivedProcEdge[eI] != -1)
                    {
                        receivedProcEdge[eI] =
                            findIndex
                            (
                                procPatch.nbrEdges(),
                                receivedProcEdge[eI]
                            );
                    }
                }
            }


            // Create a map from procEdge to receivedProcEdge
            labelList procEdgeToReceivedProcEdgeID(procEdge.size(), -1);
            forAll(procEdge, eI)
            {
                const label procEdgeID = procEdge[eI];

                // Find procEdgeID in receivedProcEdge list
                forAll(receivedProcEdge, reI)
                {
                    const label rProcEdgeID = receivedProcEdge[reI];

                    if (procEdgeID == rProcEdgeID)
                    {
                        procEdgeToReceivedProcEdgeID[eI] = reI;
                        break;
                    }
                }

                // if (procEdgeToReceivedProcEdgeID[eI] == -1)
                // {
                //     FatalErrorIn("faceCracker::detachCoupledFaces()")
                //         << "procEdgeToReceivedProcEdgeID[eI] not found" << nl
                //         << "procEdge " << procEdge << nl
                //         << "receivedProcEdge " << receivedProcEdge
                //         << abort(FatalError);
                // }
            }

            // Sync the edgeIsInternal list: an edge is internal if both
            // processors say it is internal
            forAll(procEdge, eI)
            {
                const label reI = procEdgeToReceivedProcEdgeID[eI];

                if (reI != -1)
                {
                    if (!receivedProcEdgeIsInternal[reI])
                    {
                        edgeIsInternal[eI] = false;
                    }
                }
            }
        }


        // Check global edges
        {
            const labelList sharedEdges = sharedEdgeSet.toc();

            // Note: we use a labelList instead of a boolList because there is
            // no reduce(..., orOp<boolList>(...)) function
            labelList checkSharedEdge(mesh.globalData().nGlobalEdges(), 0);

            forAll(sharedEdges, eI)
            {
                const label edgeID = sharedEdges[eI];

                const label glEdgeID = findIndex(glEdges, edgeID);

                if (glEdgeID != -1)
                {
                    checkSharedEdge[glEdgeAddr[glEdgeID]] = 1;
                }
            }

            reduce(checkSharedEdge, maxOp<labelList>());

            labelList internalSharedEdge(mesh.globalData().nGlobalEdges(), 1);

            forAll(checkSharedEdge, glEdgeI)
            {
                if (checkSharedEdge[glEdgeI])
                {
                    const label sharedEdgeIndex =
                        findIndex(glEdgeAddr, glEdgeI);

                    if (sharedEdgeIndex != -1)
                    {
                        const label curSharedEdge = glEdges[sharedEdgeIndex];

                        const labelList& curFaces = edgeFaces[curSharedEdge];

                        // Check if the edge is internal on this processor
                        forAll(curFaces, fI)
                        {
                            const label faceID = curFaces[fI];

                            if (!mesh.isInternalFace(faceID))
                            {
                                const label patchID =
                                    mesh.boundaryMesh().whichPatch(faceID);

                                if
                                (
                                   !isA<processorPolyPatch>
                                    (
                                        mesh.boundaryMesh()[patchID]
                                    )
                                )
                                {
                                    // The edge belongs to a boundary face
                                    internalSharedEdge[glEdgeI] = 0;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            reduce(internalSharedEdge, minOp<labelList>());

            forAll(sharedEdges, eI)
            {
                const label edgeID = sharedEdges[eI];

                const label glEdgeID = findIndex(glEdges, edgeID);

                if (glEdgeID != -1)
                {
                    // Update local edgeIsInternal list
                    forAll(curFaceEdges, faceToBreakEdgeI)
                    {
                        if (curFaceEdges[faceToBreakEdgeI] == edgeID)
                        {
                            if (edgeIsInternal[faceToBreakEdgeI])
                            {
                                edgeIsInternal[faceToBreakEdgeI] =
                                    internalSharedEdge[glEdgeAddr[glEdgeID]];
                                break;
                            }
                        }
                    }
                }
            }
        }


        // 2) Find added points
        //        Points that are only on a boundary edge

        labelHashSet pointsToAddSet(faceToBreak.size());

        // Initially add all points of the face, then we will remove those
        // on an internal edge
        forAll(faceToBreak, pI)
        {
            pointsToAddSet.insert(faceToBreak[pI]);
        }

        forAll(edgeIsInternal, eI)
        {
            if (edgeIsInternal[eI])
            {
                const label edgeID = curFaceEdges[eI];
                const edge& curEdge = edges[edgeID];

                // Remove points from list if found

                if (pointsToAddSet.found(curEdge.start()))
                {
                    pointsToAddSet.erase(curEdge.start());
                }

                if (pointsToAddSet.found(curEdge.end()))
                {
                    pointsToAddSet.erase(curEdge.end());
                }
            }
        }


        // Add global points to the pointsToAddSet
        labelHashSet sharedPointsToAddSet;
        {
            const labelList pointsToAddBeforeGlobalCheck = pointsToAddSet.toc();

            forAll(pointsToAddBeforeGlobalCheck, pI)
            {
                const label pointID = pointsToAddBeforeGlobalCheck[pI];
                const label glPointID = findIndex(glPoints, pointID);

                if (glPointID != -1)
                {
                    sharedPointsToAddSet.insert(pointID);
                }
            }
        }

        const labelList sharedPointsToAdd = sharedPointsToAddSet.toc();

        labelList checkSharedPoint(mesh.globalData().nGlobalPoints(), 0);

        forAll(sharedPointsToAdd, pI)
        {
            const label pointID = sharedPointsToAdd[pI];

            const label glPointID = findIndex(glPoints, pointID);

            if (glPointID != -1)
            {
                checkSharedPoint[glPointAddr[glPointID]] = 1;
            }
        }

        reduce(checkSharedPoint, maxOp<labelList>());

        forAll(checkSharedPoint, glPointI)
        {
            if (checkSharedPoint[glPointI])
            {
                const label sharedPointIndex = findIndex(glPointAddr, glPointI);

                if (sharedPointIndex != -1)
                {
                    const label pointID = glPoints[sharedPointIndex];

                    // Add point to be split
                    if (!pointsToAddSet.found(pointID))
                    {
                        pointsToAddSet.insert(pointID);
                    }
                }
            }
        }


        const labelList pointsToAdd = pointsToAddSet.toc();

        // Print points-to-add
        if (debug)
        {
            const pointField& points = mesh.points();

            Pout<< nl << "pointsToAdd: " << endl;

            forAll(pointsToAdd, pI)
            {
                Pout<< "    " << points[pointsToAdd[pI]] << endl;
            }
        }


        // 3) Check if surrounding faces need to be modified
        //        First add all point-faces of points-to-add
        //        Then, starting at the master face (face-to-break), perform a
        //        cell-face-cell walk, through faces with a point-to-add,
        //        removing faces from the faces-to-modify list

        // A processor may have a point to add, but if it does not need to
        // modify faces then we do not need to add a point
        // Try find a boundary face attached to the points-to-add; if we cannot
        // find a "master face" then we will not add points on this proc
        if (faceCellID == -1 && pointsToAdd.size())
        {
            // This proc does not have a face to break but it does need to
            // split an edge
            // Pick a boundary face to be the master face
            label masterFaceID = -1;
            forAll(pointsToAdd, pI)
            {
                const label pointID = pointsToAdd[pI];
                const labelList& curPointFaces = pointFaces[pointID];

                forAll(curPointFaces, fI)
                {
                    const label faceID = curPointFaces[fI];

                    if (!mesh.isInternalFace(faceID))
                    {
                        const label patchID =
                            mesh.boundaryMesh().whichPatch(faceID);

                        if
                        (
                            !isA<processorPolyPatch>
                            (
                                mesh.boundaryMesh()[patchID]
                            )
                        )
                        {
                            // This will be the master face
                            masterFaceID = faceID;
                            break;
                        }
                    }
                }

                if (masterFaceID != -1)
                {
                    break;
                }
            }

            if (masterFaceID != -1)
            {
                faceCellID = faceOwn[masterFaceID];
            }
        }

        // First add all pointFaces of pointsToAdd to the list of faces to check
        labelHashSet facesToModifySet(20);

        if (faceCellID != -1)
        {
            forAll(pointsToAdd, pI)
            {
                const labelList& curPointFaces = pointFaces[pointsToAdd[pI]];

                forAll(curPointFaces, pfI)
                {
                    const label pointFaceID = curPointFaces[pfI];

                    if (!facesToModifySet.found(pointFaceID))
                    {
                        facesToModifySet.insert(pointFaceID);
                    }
                }
            }

            // Starting at the faceCell, check all faces of the cell; if the
            // face contains a point, we will remove it from facesToModifySet
            // and add the neighbour cell to be checked next. We will continue
            // this cell-face-cell walk until there are no more cells to check.

            if (faceCellID != -1)
            {
                SLList<label> cellsToCheck;
                //const label faceCellID = faceOwn[faceToBreakID];
                cellsToCheck.append(faceCellID);
                labelHashSet checkedFaces(20);

                const cellList& cells = mesh.cells();

                // Cell-face-cell walk
                do
                {
                    const label cellID = cellsToCheck.first();
                    const labelList& curCellFaces = cells[cellID];

                    forAll(curCellFaces, fI)
                    {
                        const label faceID = curCellFaces[fI];

                        if (!checkedFaces.found(faceID))
                        {
                            checkedFaces.insert(faceID);

                            // Remove from the facesToModifySet if found
                            if (facesToModifySet.found(faceID))
                            {
                                facesToModifySet.erase(faceID);

                                // Add neighbour cell to cellsToCheck
                                if (mesh.isInternalFace(faceID))
                                {
                                    label neiCellID = faceNei[faceID];
                                    if (neiCellID == cellID)
                                    {
                                        neiCellID = faceOwn[faceID];
                                    }

                                    cellsToCheck.append(neiCellID);
                                }
                            }
                        }
                    }

                    // Remove current cell from cellsToCheck
                    cellsToCheck.removeHead();
                }
                while (cellsToCheck.size());
            }
        }

        if (facesToModifySet.size())
        {
            const labelList facesToModify = facesToModifySet.toc();

            // Print faces to be modified
            if (debug)
            {
                Pout<< nl << "facesToModify: " << endl;

                forAll(facesToModify, fI)
                {
                    Pout<< "Modify face "
                        << mesh.faceCentres()[facesToModify[fI]]
                        << endl;
                }
            }


            // 4) If faces were modified
            //        Add the points-to-add

            labelList addedPoints(pointsToAdd.size(), -1);
            const pointField& points = mesh.points();

            forAll(pointsToAdd, pI)
            {
                const label pointID = pointsToAdd[pI];

                addedPoints[pI] =
                    ref.setAction
                    (
                        polyAddPoint
                        (
                            points[pointID],           // point
                            pointID,                   // master point
                            -1,                        // zone ID
                            true                       // supports a cell
                        )
                    );
            }


            // 5) Modify surrounding faces

            forAll(facesToModify, fI)
            {
                const label faceID = facesToModify[fI];

                const face& oldFace = faces[faceID];
                face newFace = oldFace;

                forAll (oldFace, pI)
                {
                    const label oldPointID = oldFace[pI];

                    // Check if the point has been split
                    forAll(pointsToAdd, ptoaI)
                    {
                        const label pointToAddID = pointsToAdd[ptoaI];

                        if (pointToAddID == oldPointID)
                        {
                            newFace[pI] = addedPoints[ptoaI];
                            break;
                        }
                    }
                }

                if (mesh.isInternalFace(faceID))
                {
                    ref.setAction
                    (
                        polyModifyFace
                        (
                            newFace,                    // face
                            faceID,                     // master face
                            faceOwn[faceID],            // owner
                            faceNei[faceID],            // neighbour
                            false,                      // flip flux
                            -1,                         // patch for face
                            false,                      // remove from zone
                            -1,                         // zone for face
                            false                       // face zone flip
                        )
                    );
                }
                else
                {
                    ref.setAction
                    (
                        polyModifyFace
                        (
                            newFace,                     // face
                            faceID,                      // master face
                            faceOwn[faceID],             // owner
                            -1,                          // neighbour
                            false,                       // flip flux
                            mesh.boundaryMesh().whichPatch(faceID), // patch
                            false,                       // remove from zone
                            -1,                          // zone for face
                            false                        // face zone flip
                        )
                    );
                }
            }
        }


        // 6) Move the coupled face to the crack patch

        if (faceToBreakID != -1)
        {
            ref.setAction
            (
                polyModifyFace
                (
                    faceToBreak,                 // face
                    faceToBreakID,               // master face
                    faceCellID,                  // owner
                    -1,                          // neighbour
                    false,                       // flip flux
                    crackPatchID_.index(),       // patch
                    false,                       // remove from zone
                    -1,                          // zone for face
                    false                        // face zone flip
                )
            );
        }

        // Clear the coupled-faces-to-break list
        coupledFacesToBreak_.setSize(0);
    }
}


// ************************************************************************* //
