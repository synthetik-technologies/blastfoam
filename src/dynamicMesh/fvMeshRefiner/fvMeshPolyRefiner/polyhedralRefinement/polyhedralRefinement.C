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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "polyhedralRefinement.H"
#include "polyTopoChange.H"
#include "polyAddCell.H"
#include "polyAddPoint.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "meshTools.H"
#include "foamMeshTools.H"
#include "syncTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyhedralRefinement, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::polyhedralRefinement::getAnchorLevel
(
    const label faceI
) const
{
    const face& f = mesh_.faces()[faceI];

    if (f.size() <= 3)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        const label& ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) >= 3)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel + 1) >= 3)
        {
            return ownLevel + 1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::polyhedralRefinement::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label cellI,

    polyTopoChange& ref
) const
{
    // Find in every face the cellLevel + 1 points (from edge subdivision)
    // and the anchor points

    // Get current cell and its level
    const cell& curCell = mesh_.cells()[cellI];
    const label& cLevel = cellLevel_[cellI];

    // Get mesh faces and face edges
    const faceList& meshFaces = mesh_.faces();
    const labelListList& meshFaceEdges = mesh_.faceEdges();

    // Get mesh cell points
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // Map from edge mid to anchor points
    Map<edge> midPointToAnchors(24);
    // Map from edge mid to face mids
    Map<edge> midPointToFaceMids(24);

    // Running count of number of internal faces added so far
    label nFacesAdded = 0;

    // Loop through faces of the cell
    forAll(curCell, i)
    {
        // Get face index
        const label& faceI = curCell[i];

        // Get current face and its edges
        const face& f = meshFaces[faceI];
        const labelList& fEdges = meshFaceEdges[faceI];

        // We are on the cellI side of face f. The face will have 1 or n
        // cLevel points (where n is the number of points/edges of a face)
        // and lots of higher numbered ones

        // Index of face mid point
        label faceMidPointI = -1;

        // Get number of anchors for the face
        const label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel + 1 and cLevel + 2 points

            // Find the one anchor
            label anchorFp = -1;

            // Loop through face points
            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    // Point level is smaller than cLevel + 1, this is the
                    // anchor point
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel + 1 point
            label edgeMid = findLevel(f, f.fcIndex(anchorFp), true, cLevel + 1);
            label faceMid = findLevel(f, f.fcIndex(edgeMid), true, cLevel + 1);

            // Set face mid point index
            faceMidPointI = f[faceMid];
        }
        else
        {
            // There is no face middle yet but the face will be split. Set face
            // mid point index
            faceMidPointI = faceMidPoint[faceI];
        }


        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointInfo will collect
        // all the info and combine it all
        forAll(f, fp0)
        {
            // Get point index
            const label& point0 = f[fp0];

            if (pointLevel_[point0] <= cLevel)
            {
                // This is anchor point

                // Walk forward to cLevel + 1 or edgeMidPoint of this level
                label edgeMidPointI = -1;

                const label fp1 = f.fcIndex(fp0);

                if (pointLevel_[f[fp1]] <= cLevel)
                {
                    // Another anchor: edge will be split
                    const label& edgeI = fEdges[fp0];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    // Sanity check
                    if (edgeMidPointI == -1)
                    {
                        const labelList& cPoints = meshCellPoints[cellI];

                        FatalErrorInFunction
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search forward in face to clevel + 1
                    const label edgeMid = findLevel(f, fp1, true, cLevel + 1);

                    edgeMidPointI = f[edgeMid];
                }

                label newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    true,                   // mid point after anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    ref
                );

                if (newFaceI != -1)
                {
                    ++nFacesAdded;
                }


                // Now walk backward

                label fpMin1 = f.rcIndex(fp0);

                if (pointLevel_[f[fpMin1]] <= cLevel)
                {
                    // Another anchor: edge will be split
                    const label& edgeI = fEdges[fpMin1];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    // Sanity check
                    if (edgeMidPointI == -1)
                    {
                        const labelList& cPoints = meshCellPoints[cellI];

                        FatalErrorInFunction
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search back in face to clevel + 1
                    const label edgeMid = findLevel(f, fpMin1, false, cLevel + 1);

                    edgeMidPointI = f[edgeMid];
                }

                newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    false,                  // mid point before anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    ref
                );

                if (newFaceI != -1)
                {
                    ++nFacesAdded;
                }
            } // End for this anchor point
        } // End for all face points
    } // End for all cell faces
}


Foam::label Foam::polyhedralRefinement::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label cellI,
    const label faceI,
    const label pointI
) const
{
    // Check whether pointI is an anchor of cellI. If it is not, check whether
    // any other point on the face is an anchor for the cell

    if (cellAnchorPoints[cellI].size() > 0)
    {
        const label index = findIndex(cellAnchorPoints[cellI], pointI);

        if (index != -1)
        {
            return cellAddedCells[cellI][index];
        }


        // pointI is not an anchor for the cell. Maybe we already refined the
        // face so check all the face vertices
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            const label index = findIndex(cellAnchorPoints[cellI], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[cellI][index];
            }
        }

        // Problem: the point does not seem to be an anchor for cell

        // Pick up points of the cell
        const labelList& cPoints = mesh_.cellPoints()[cellI];

        Perr<< "cell: " << cellI << ", points: " << endl;
        forAll(cPoints, i)
        {
            const label pointI = cPoints[i];

            Perr<< "    " << pointI << " coord: " << mesh_.points()[pointI]
                << nl;
        }

        Perr<< "cell: " << cellI << " anchorPoints: " << cellAnchorPoints[cellI]
            << endl;

        FatalErrorInFunction
            << "Could not find point " << pointI
            << " in the anchorPoints for cell " << cellI << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return cellI;
    }
}


void Foam::polyhedralRefinement::setNewFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label faceI,
    const label pointI,

    label& own,
    label& nei
) const
{
    // Get anchor cell for this anchor point on owner side
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[faceI],
        faceI,
        pointI
    );

    if (mesh_.isInternalFace(faceI))
    {
        // Get anchor cell for this anchor point on neighbour side
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[faceI],
            faceI,
            pointI
        );
    }
    else
    {
        // Boundary face: set neighbour to -1
        nei = -1;
    }
}


Foam::label Foam::polyhedralRefinement::findLevel
(
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointI = f[fp];

        if (pointLevel_[pointI] < wantedLevel)
        {
            FatalErrorInFunction
                << "face:" << f
                << " level:" << IndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointI] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    FatalErrorInFunction
        << "face:" << f
        << " level:" << IndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


Foam::label Foam::polyhedralRefinement::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label cellI,
    const label faceI,
    const bool faceOrder,
    const label edgeMidPointI,
    const label anchorPointI,
    const label faceMidPointI,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& ref
) const
{
    // A single internal face is added per edge inbetween anchor points,
    // i.e. one face per midPoint between anchor points. The information is
    // stored on the midPoint and if we have enough information (finished
    // collecting two anchors and two face mid points), we add the face.
    // Note that this member function can get called anywhere from
    // two times (two unrefined faces) to four times (two refined faces) so
    // the first call that adds the information creates the face


    // See if need to store anchors
    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointI);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointI, edge(anchorPointI, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointI);

    if (faceMidFnd == midPointToFaceMids.end())
    {
        midPointToFaceMids.insert(edgeMidPointI, edge(faceMidPointI, -1));
    }
    else
    {
        edge& e = faceMidFnd();

        if (faceMidPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = faceMidPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Check if this call of storeMidPointInfo is the one that completed all
    // the nessecary information

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        const edge& anchors = midPointToAnchors[edgeMidPointI];
        const edge& faceMids = midPointToFaceMids[edgeMidPointI];

        label otherFaceMidPointI = faceMids.otherVertex(faceMidPointI);

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cell mid and face mids

        DynamicList<label> newFaceVerts(4);
        if (faceOrder == (mesh_.faceOwner()[faceI] == cellI))
        {
            newFaceVerts.append(faceMidPointI);

            // Check and insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointI,  // edge between faceMid and
                edgeMidPointI,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                otherFaceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }
        else
        {
            newFaceVerts.append(otherFaceMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointI,
                edgeMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                faceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }

        face newFace;
        newFace.transfer(newFaceVerts.shrink());
        newFaceVerts.clear();

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchorPointI
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchors.otherVertex(anchorPointI)
        );

        // Get mesh points
        const pointField& meshPoints = mesh_.points();

        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = meshPoints[anchorPointI];
            neiPt = meshPoints[anchors.otherVertex(anchorPointI)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace = newFace.reverseFace();

            ownPt = meshPoints[anchors.otherVertex(anchorPointI)];
            neiPt = meshPoints[anchorPointI];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = meshPoints[anchorPointI];
                neiPt = meshPoints[anchors.otherVertex(anchorPointI)];
            }
            else
            {
                ownPt = meshPoints[anchors.otherVertex(anchorPointI)];
                neiPt = meshPoints[anchorPointI];
            }

            checkInternalOrientation
            (
                ref,
                cellI,
                faceI,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            ref,
            faceI,
            anchorPointI,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}


void Foam::polyhedralRefinement::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    // Get number of points
    const label nPoints = mesh_.nPoints();

    if (p0 < nPoints && p1 < nPoints)
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::polyhedralRefinement::setRefinement
(
    polyTopoChange& meshMod,
    const labelList& cellsToRefine
) const
{
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, celli)
    {
        newCellLevel.append(cellLevel_[celli]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointi)
    {
        newPointLevel.append(pointLevel_[pointi]);
    }

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Allocating " << cellsToRefine.size() << " cell midpoints."
            << endl;
    }


    // PART 1: Mark cells for refinement and add points at their cell centres

    // Get necessary mesh data
    const faceList& meshFaces = mesh_.faces();
    const cellList& meshCells = mesh_.cells();
    const vectorField& meshCellCentres = mesh_.cellCentres();

    // Mid point for refined cell (points at cell centres):
    // Not refined = -1
    // Shall be refined > -1 (label of added mid point)
    labelList cellMidPoint(mesh_.nCells(), -1);

    // Loop through cells to refine
    forAll(cellsToRefine, i)
    {
        // Get cell idnex
        const label& cellI = cellsToRefine[i];
        label anchorPointi = mesh_.faces()[mesh_.cells()[cellI][0]][0];

        cellMidPoint[cellI] = meshMod.setAction
        (
            polyAddPoint
            (
                meshCellCentres[cellI], // Point to add (cell centre)
                anchorPointi,           // Appended point: no master ID
                -1,                     // Zone for point
                true                    // Supports a cell
            )
        );
        newPointLevel(cellMidPoint[cellI]) = cellLevel_[cellI] + 1;
    }

    // Write out split cells as a cell set for debug
    if (debug)
    {
        // Note: cellSet is actually a hash table of labels
        cellSet splitCells(mesh_, "splitCells", cellsToRefine.size());

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] > -1)
            {
                // Cell is marked for refinement, insert into cellSet
                splitCells.insert(cellI);
            }
        }

        Pout<< FUNCTION_NAME << nl
            << "Writing " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }


    // PART 2: Mark edges for refinement and add points to edge centres

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Allocating edge midpoints."
            << endl;
    }

    // Refined edges are defined by having both their point levels lower than
    // the cell level, i.e. if any cell that gets split uses this edge, the edge
    // needs to be split as well

    // Get necessary mesh data
    const labelListList& meshCellEdges = mesh_.cellEdges();
    const edgeList& meshEdges = mesh_.edges();

    // Mid points for refined edge:
    // No need to split edge = -1
    // Label of introduced mid point > -1
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over all cells instead of all edges
    forAll(cellMidPoint, cellI)
    {
        // Point is marked for refinement, proceed to look at the edges
        if (cellMidPoint[cellI] > -1)
        {
            // Get edges of this cell
            const labelList& cEdges = meshCellEdges[cellI];

            forAll(cEdges, i)
            {
                // Get edge index and edge
                const label& edgeI = cEdges[i];
                const edge& e = meshEdges[edgeI];

                if
                (
                    pointLevel_[e[0]] <= cellLevel_[cellI]
                 && pointLevel_[e[1]] <= cellLevel_[cellI]
                )
                {
                    // Point levels of both edge points are <= cell level, mark
                    // the edge for splitting
                    edgeMidPoint[edgeI] = 12345;
                }
            }
        }
    }

    // Synchronize edgeMidPoint across coupled patches. Note: use max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
//         false               // no separation
    );

    // Introduce edge points

    // Get necessary mesh data
    const pointField& meshPoints = mesh_.points();

    // Memory management
    {
        // Phase 1: calculate midpoints and sync. This is necessary if we don't
        // use binary format for writing and we slowly get differences.

        // Allocate storage for edge points
        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                // Edge marked to be split. Get edge centre
                edgeMids[edgeI] = meshEdges[edgeI].centre(meshPoints);
            }
        }

        // Sync across processor boundaries
        syncTools::syncEdgeList
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
//             true                           // apply separation
        );

        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                const edge& e = mesh_.edges()[edgeI];
                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI], // Point
                        e[0],            // Appended point, no master ID
                        -1,              // Zone for point
                        true             // Supports a cell
                    )
                );
                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    ) + 1;
            }
        }
    } // End memory management for syncing/adding edge points

    // Write out edge mid points for split edges for debugging
    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                // Get edge and write its cell centre
                const edge& e = meshEdges[edgeI];
                meshTools::writeOBJ(str, e.centre(meshPoints));
            }
        }

        Pout<< FUNCTION_NAME << nl
            << "Writing centres of edges to split to file " << str.name()
            << endl;
    }


    // PART 3: Calculate face level (after selected cells splitting)

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed at least 3 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); ++faceI)
    {
        faceAnchorLevel[faceI] = getAnchorLevel(faceI);
    }

    // Mid points for refined face (points at face centres):
    // Not refined = -1
    // Shall be refined > -1 (label of added mid point)
    labelList faceMidPoint(mesh_.nFaces(), -1);

    // Get necessary mesh data
    const labelList& meshFaceOwner = mesh_.faceOwner();
    const labelList& meshFaceNeighbour = mesh_.faceNeighbour();

    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Internal faces: look at cells on both sides. Uniquely determined since
    // the face itself is guaranteed to be same level as most refined neighbour
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Note: no need to check whether the face has valid anchor level since
        // all faces can be split
        const label& own = meshFaceOwner[faceI];
        const label& ownLevel = cellLevel_[own];
        const label newOwnLevel = ownLevel + (cellMidPoint[own] > -1 ? 1 : 0);

        const label& nei = meshFaceNeighbour[faceI];
        const label& neiLevel = cellLevel_[nei];
        const label newNeiLevel = neiLevel + (cellMidPoint[nei] > -1 ? 1 : 0);

        if
        (
            newOwnLevel > faceAnchorLevel[faceI]
         || newNeiLevel > faceAnchorLevel[faceI]
        )
        {
            // New level is higher than the face anchor level, mark for
            // splitting
            faceMidPoint[faceI] = 12345;
        }
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    // refinining and subsetting)

    // Memory management
    {
        // Create list for swapping boundary data
        labelList newNeiLevel(nFaces - nInternalFaces);

        forAll(newNeiLevel, i)
        {
            const label& own = meshFaceOwner[i + nInternalFaces];
            const label& ownLevel = cellLevel_[own];
            const label newOwnLevel =
                ownLevel + (cellMidPoint[own] > -1 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap the list which now contains data from the other side
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);
//         syncTools::swapBoundaryFaceList(mesh_, newNeiLevel, false);

        forAll(newNeiLevel, i)
        {
            // Get face index
            const label faceI = i + nInternalFaces;

            // Note: no need to check whether the face has valid anchor level
            // since all faces can be split
            const label& own = meshFaceOwner[faceI];
            const label& ownLevel = cellLevel_[own];
            const label newOwnLevel =
                ownLevel + (cellMidPoint[own] > -1 ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[faceI]
             || newNeiLevel[i] > faceAnchorLevel[faceI]
            )
            {
                // New level is higher than the face anchor level, mark for
                // splitting
                faceMidPoint[faceI] = 12345;
            }
        }
    } // End memory management for syncing owner/neighbour face levels

    // Note: synronisation of faceMidPoints across coupled patches is not
    // necessary since we have exchanged the neighbour data above using
    // swapBoundaryFaceList, thus the faceMidPoint has to be the same on both
    // sides. VV, 5/Jan/2018.
//    // Synchronize faceMidPoint across coupled patches
//    syncTools::syncFaceList
//    (
//        mesh_,
//        faceMidPoint,
//        maxEqOp<label>(),
//        false
//    );


    // Introduce face points

    // Get face centres
    const vectorField& meshFaceCentres = mesh_.faceCentres();

    // Memory management
    {
        // Phase 1: determine mid points and sync. Note: the same procedure has
        // been used for syncing edge mid points

        // Allocate storage for boundary face points
        pointField bFaceMids
        (
            nFaces - nInternalFaces,
            point(-GREAT, -GREAT, -GREAT)
        );

        // Loop through boundary face mids
        forAll(bFaceMids, i)
        {
            // Get face index
            const label faceI = i + nInternalFaces;

            if (faceMidPoint[faceI] > -1)
            {
                // This is a valid face mid, get the face centre
                bFaceMids[i] = meshFaceCentres[faceI];
            }
        }

        // Sync across coupled boundaries. Note: uses maximum of the components
        // of the vector. This is completely arbitrary but it doesn't matter as
        // long as we have same points on each side. VV, 5/Jan/2018.
        syncTools::syncBoundaryFaceList
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
//             true               // apply separation
        );

        // Loop through faces
        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                const face& f = mesh_.faces()[faceI];
                // Face marked to be split. Add the point at face centre and
                // replace faceMidPoint with actual point label

                faceMidPoint[faceI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            faceI < nInternalFaces
                          ? meshFaceCentres[faceI]
                          : bFaceMids[faceI - nInternalFaces]
                        ),    // Point
                        f[0], // Appended point, no master ID
                        -1,   // Zone for point
                        true  // Supports a cell
                    )
                );
                newPointLevel(faceMidPoint[faceI]) = faceAnchorLevel[faceI]+1;
            }
        }
    } // End memory management for syncing boundary data and adding face mids

    // Write out split faces as a face set for debugging
    if (debug)
    {
        // Create a faceSet with 3*cell sizes to prevent excessive resizing
        faceSet splitFaces(mesh_, "splitFaces", 3*cellsToRefine.size());

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                splitFaces.insert(faceI);
            }
        }

        Pout<< FUNCTION_NAME << nl
            << "Writing " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath()
            << endl;

        splitFaces.write();
    }


    // Now we have all the information we need to perform the refinement and we
    // no longer need to refer to cellsToRefine_. The information is:
    // - cellMidPoint >= 0 : cell needs to be split
    // - faceMidPoint >= 0 : face needs to be split
    // - edgeMidPoint >= 0 : edge needs to be split


    // PART 4: Get corner and anchor points for all cells

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Finding cell anchorPoints" << endl;
    }

    // Get anchor points for each cell: points that have the same or lower
    // refinement level as the cell
    List<DynamicList<label>> cellAnchorPointsDynamic(mesh_.nCells());

    // Loop through all cells
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            // This cell shall be cut. Set capacity to 16 to prevent excessive
            // resizing
            cellAnchorPointsDynamic[cellI].setCapacity(16);
        }
    }

    // Get necessary mesh data
    const labelListList& meshPointCells = mesh_.pointCells();

    // Loop through all points
    forAll(pointLevel_, pointI)
    {
        // Get point cells
        const labelList& pCells = meshPointCells[pointI];

        // Loop through all cells sharing this point
        forAll(pCells, pCellI)
        {
            // Get current cell index
            const label& cellI = pCells[pCellI];

            if
            (
                cellMidPoint[cellI] > -1
             && pointLevel_[pointI] <= cellLevel_[cellI]
            )
            {
                // This point cells is marked for refinement and its point level
                // is smaller or equal to cell level, append the point
                cellAnchorPointsDynamic[cellI].append(pointI);
            }
        }
    }

    // Loop through all cells and check whether at least 4 anchor points
    // have been found (minimum requirement for a tet cell)

    // Get cell points for error output
    const labelListList& meshCellPoints = mesh_.cellPoints();

    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            // Cell selected for refinement
            if (cellAnchorPointsDynamic[cellI].size() < 4)
            {
                // Cell has less than 4 anchor points. Issue an error and report
                // cell points
                const labelList& cPoints = meshCellPoints[cellI];

                FatalErrorInFunction
                   << "Cell " << cellI
                    << " of level " << cellLevel_[cellI]
                    << " does not seem to have enough points of "
                    << " lower level" << endl
                    << "cellPoints:" << cPoints << endl
                    << "pointLevels:"
                    << IndirectList<label>(pointLevel_, cPoints)() << endl
                    << abort(FatalError);
            }
        }
    }

    // Collect cellAnchorPoints into a List<labelList> instead of
    // List<dynamicList>
    labelListList cellAnchorPoints(mesh_.nCells());

    forAll(cellAnchorPointsDynamic, cellI)
    {
        // Tranfer the dynamic list for each cell into an ordinary list
        cellAnchorPoints[cellI].transfer(cellAnchorPointsDynamic[cellI]);
    }

    // PART 5: Add the cells

    if (debug)
    {
        Pout<< "polyhedralRefinement::setRefinementInstruction(...)" << nl
            << " Adding cells."
            << endl;
    }

    // We should have exactly n new cells per each split cell, where n is the
    // number of anchor points in a cell
    labelListList cellAddedCells(mesh_.nCells());

    // Get cell zone mesh
    const meshCellZones& cellZones = mesh_.cellZones();

    forAll(cellAnchorPoints, cellI)
    {
        // Check whether this is a split cell
        if (cellMidPoint[cellI] > -1)
        {
            // Get cell anchors
            const labelList& cAnchors = cellAnchorPoints[cellI];

            // Set the total number of added cells to number of anchors
            labelList& cAdded = cellAddedCells[cellI];
            cAdded.setSize(cAnchors.size());

            // Original cell has index 0
            cAdded[0] = cellI;

            // Update cell level
            newCellLevel[cellI] = cellLevel_[cellI]+1;

            // Add other cells
            for (label i = 1; i < cAdded.size(); ++i)
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                         // Master point
                        -1,                         // Master edge
                        -1,                         // Master face
                        cellI,                      // Master cell
                        cellZones.whichZone(cellI)  // Zone for cell
                    )
                );
                newCellLevel(cAdded[i]) = cellLevel_[cellI]+1;
            }
        }
    }


    // PART 6: Adding faces

    // 6.1. Existing faces that get split (into n faces where n is the number of
    //      points or edges)
    // 6.2. Existing faces that do not get split but only edges get split
    // 6.3. Existing faces that do not get split but get new owner/neighbour
    // 6.4. New internal faces inside split cells.

    if (debug)
    {
        Pout<< "polyhedralRefinement::setRefinementInstruction(...)" << nl
            << " Marking faces to be handled"
            << endl;
    }

    // Get all faces to split:
    // a) All faces of a cell being split
    // b) All faces that are being split
    // c) Both faces of an edge that is being split
    boolList facesToSplit(mesh_.nFaces(), false);

    // Get edge faces
    const labelListList& meshEdgeFaces = mesh_.edgeFaces();

    // a) All faces of a cell that is being split
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            const cell& cFaces = meshCells[cellI];

            forAll(cFaces, i)
            {
                facesToSplit[cFaces[i]] = true;
            }
        }
    }

    // b) All faces that are being split
    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] > -1)
        {
            facesToSplit[faceI] = true;
        }
    }

    // c) Both faces of an edge that are being split
    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] > -1)
        {
            const labelList& eFaces = meshEdgeFaces[edgeI];

            forAll(eFaces, i)
            {
                facesToSplit[eFaces[i]] = true;
            }
        }
    }

    // Note: after splitting a certain face during parts 6.1. to 6.4.,
    // facesToSplit for that face will be set back to 0, i.e. marked as finished


    // PART 6.1. Add/modify faces for each face being split

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Splitting faces..." << endl;
    }

    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] > -1 && facesToSplit[faceI])
        {
            // Face has not been split.
            // Note: although facesToSplit can't be different than 1 here and
            // the second check can be ommitted, it is left for clarity

            // Get the face
            const face& f = meshFaces[faceI];

            // Flag to control whether the original faceI has been used
            // Note: original face gets modified, other n - 1 faces are added,
            // where n is the number of points/edges of a face
            bool modifiedFace = false;

            // Get anchor level for the face
            const label& anchorLevel = faceAnchorLevel[faceI];

            // New face always has four points/edges for arbitrary polygon
            face newFace(4);

            // Loop through all points of original face
            forAll(f, fp)
            {
                const label& pointI = f[fp];

                if (pointLevel_[pointI] <= anchorLevel)
                {
                    // This point is anchor, start collecting face

                    // Create dynamic list (because of append) for face vertices
                    // and append the first (anchor) point
                    DynamicList<label> faceVerts(4);
                    faceVerts.append(pointI);

                    // Walk forward to mid point.
                    // - if next is +2 midpoint is +1
                    // - if next is +1 it is midpoint
                    // - if next is +0 there has to be edgeMidPoint

                    // Appends all points from this point to face mid point
                    walkFaceToMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    // Append face mid point
                    faceVerts.append(faceMidPoint[faceI]);

                    // Append all points from face mid point to starting point
                    walkFaceFromMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    // Transfer dynamic list to a face (ordinary list)
                    newFace.transfer(faceVerts);
                    faceVerts.clear();

                    // Set new owner/neighbour indices based on split cells
                    label own, nei;
                    setNewFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        pointI,           // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        // Get mesh cell centres
                        const vectorField& meshCellCentres =
                            mesh_.cellCentres();

                        if (mesh_.isInternalFace(faceI))
                        {
                            const label oldOwn = meshFaceOwner[faceI];
                            const label oldNei = meshFaceNeighbour[faceI];

                            // Print info only with deep debug level
                            if (debug > 1)
                            {
                                Pout<< "Split internal face: " << faceI
                                    << ", verts: " << f
                                    << ", into quad: " << newFace << nl
                                    << "owner: " << oldOwn
                                    << ", neighbour: " << oldNei << endl;
                            }

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                meshCellCentres[oldOwn],
                                meshCellCentres[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            const label oldOwn = meshFaceOwner[faceI];

                            // Print info only with deep debug level
                            if (debug > 1)
                            {
                                Pout<< "Split boundary face: " << faceI
                                    << ", verts: " << f
                                    << ", into quad: " << newFace << nl
                                    << "owner: " << oldOwn << endl;
                            }

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                meshCellCentres[oldOwn],
                                meshFaceCentres[faceI],
                                newFace
                            );
                        }
                    } // End debug


                    if (!modifiedFace)
                    {
                        // Modify first face
                        modifiedFace = true;
                        modifyFace(meshMod, faceI, newFace, own, nei);
                    }
                    else
                    {
                        // Add additional faces
                        addFace(meshMod, faceI, newFace, own, nei);
                    }
                } // End point anchor check
            } // End for all points

            // Mark face as handled
            facesToSplit[faceI] = false;
        }
    }


    // PART 6.2. Modify faces that do not get split but have edges that are
    // being split

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Modifying faces with split edges..."
            << endl;
    }

    // Get face edges
    const labelListList& meshFaceEdges = mesh_.faceEdges();

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] > -1)
        {
            // This is an edge that is going to be split, get edge faces
            const labelList& eFaces = meshEdgeFaces[edgeI];

            // Loop through all faces of an edge
            forAll(eFaces, i)
            {
                // Get face index
                const label faceI = eFaces[i];

                // Check whether this is not a face that has been marked for
                // splitting and that the face has not been handled yet. The
                // second check is necessary since we go through edge faces
                // instead of just faces
                if (faceMidPoint[faceI] < 0 && facesToSplit[faceI])
                {
                    // This is unsplit face that has not been handled

                    // Get face and face edges
                    const face& f = meshFaces[faceI];
                    const labelList& fEdges = meshFaceEdges[faceI];

                    // Create a dynamic list containing new face vertices
                    DynamicList<label> newFaceVerts(f.size());

                    // Append all original points and all edge mid points
                    forAll(f, fp)
                    {
                        newFaceVerts.append(f[fp]);

                        const label edgeI = fEdges[fp];

                        if (edgeMidPoint[edgeI] > -1)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    face newFace;
                    newFace.transfer(newFaceVerts.shrink());


                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    const label anchorFp = findMinLevel(f);

                    label own, nei;
                    setNewFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        f[anchorFp],      // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        // Get mesh cell centres
                        const vectorField& meshCellCentres =
                            mesh_.cellCentres();

                        if (mesh_.isInternalFace(faceI))
                        {
                            const label oldOwn = meshFaceOwner[faceI];
                            const label oldNei = meshFaceNeighbour[faceI];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                meshCellCentres[oldOwn],
                                meshCellCentres[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            const label oldOwn = meshFaceOwner[faceI];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                meshCellCentres[oldOwn],
                                meshFaceCentres[faceI],
                                newFace
                            );
                        }
                    } // End debug

                    // Modify the face
                    modifyFace(meshMod, faceI, newFace, own, nei);

                    // Mark face as handled
                    facesToSplit[faceI] = false;

                }// End if unsplit, unhandled face
            } // End for all edge faces
        } // End if edge has been cut
    } // End for all edges


    // PART 6.3.: Modify faces that do not get split but whose owner/neighbour
    // change due to splitting

    if (debug)
    {
        Pout<< FUNCTION_NAME
            << " Changing owner/neighbour for otherwise unaffected faces..."
            << endl;
    }

    forAll(facesToSplit, faceI)
    {
        // All remaining unnaffected faces are the ones whose owner/neighbour
        // changed
        if (facesToSplit[faceI])
        {
            // Get the face
            const face& f = meshFaces[faceI];

            // The point with the lowest level should be an anchor point of the
            // neighbouring cells
            label anchorFp = findMaxLevel(f);

            label own, nei;
            setNewFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                faceI,
                f[anchorFp],      // Anchor point

                own,
                nei
            );

            // Modify the face, changing owner and neighbour
            modifyFace(meshMod, faceI, f, own, nei);

            // Mark face as handled
            facesToSplit[faceI] = false;
        }
    }


    // PART 6.4. Add new internal faces inside split cells

    // We have to find the splitting points between the anchor points. But the
    // edges between the anchor points might have been split (into two, three or
    // four edges)

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << " Adding new internal faces for split cells..."
            << endl;
    }

    // Loop through cells
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            // Cell has been split, create internal faces
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                cellI,
                meshMod
            );
        }
    }

    // Debug: check minimum point index of added points, needs to be equal to
    // number of points in the original mesh
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] > -1)
            {
                minPointI = min(minPointI, cellMidPoint[cellI]);
                maxPointI = max(maxPointI, cellMidPoint[cellI]);
            }
        }
        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                minPointI = min(minPointI, faceMidPoint[faceI]);
                maxPointI = max(maxPointI, faceMidPoint[faceI]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                minPointI = min(minPointI, edgeMidPoint[edgeI]);
                maxPointI = max(maxPointI, edgeMidPoint[edgeI]);
            }
        }

        if (minPointI != labelMax && minPointI != mesh_.nPoints())
        {
            FatalErrorInFunction
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointI: " << minPointI
                << " maxPointI: " << maxPointI
                << abort(FatalError);
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);
}


void Foam::polyhedralRefinement::setUnrefinement
(
    polyTopoChange& meshMod,
    const labelList& splitPointsToUnrefine
) const
{
    // Get point cells necessary for debug and face removal
    const labelListList& meshPointCells = mesh_.pointCells();

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Checking validity of cellLevel before setting unrefinement."
            << endl;

        forAll(cellLevel_, cellI)
        {
            if (cellLevel_[cellI] < 0)
            {
                FatalErrorInFunction
                    << "Illegal cell level " << cellLevel_[cellI]
                    << " for cell " << cellI
                    << abort(FatalError);
            }
        }

        // Write split points into a point set
        pointSet pSet
        (
            mesh_,
            "splitPoints",
            labelHashSet(splitPointsToUnrefine)
        );
        pSet.write();

        // Write split point cells into a cell set
        cellSet cSet
        (
            mesh_,
            "splitPointCells",
            splitPointsToUnrefine.size()
        );

        forAll(splitPointsToUnrefine, i)
        {
            // Get point cells and insert them into cell set
            const labelList& pCells = meshPointCells[splitPointsToUnrefine[i]];

            forAll(pCells, j)
            {
                cSet.insert(pCells[j]);
            }
        }
        cSet.write();

        Pout<< FUNCTION_NAME << nl
            << "Writing " << pSet.size()
            << " points and "
            << cSet.size() << " cells for unrefinement to" << nl
            << "pointSet " << pSet.objectPath() << nl
            << "cellSet " << cSet.objectPath()
            << endl;
    }

    // Create lists needed by face remover
    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    // Memory management
    {
        // Collect split faces in the hash set, guess size to prevent excessive
        // resizing
        labelHashSet splitFaces(12*splitPointsToUnrefine.size());

        // Get point faces
        const labelListList& meshPointFaces = mesh_.pointFaces();

        forAll(splitPointsToUnrefine, i)
        {
            // Loop through all faces of this point and insert face index
            const labelList& pFaces = meshPointFaces[splitPointsToUnrefine[i]];

            forAll(pFaces, j)
            {
                splitFaces.insert(pFaces[j]);
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // Pierced faces

            cellRegion,         // Region merged into (-1 for no region)
            cellRegionMaster,   // Master cell for region
            facesToRemove       // List of faces to be removed
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorInFunction
                << "Either the initial set of split points to unrefine does not"
                << " seem to be consistent or there are no mid points of"
                << " refined cells."
                << abort(FatalError);
        }
    }

    // Find point region master for every cell region.  This is the central point
    // from which the coarse cell will be made
    // The property of the point region master is that all cells that touch it
    // have the same cell region index
    // HJ, 6/Sep/2019
    labelList pointRegionMaster(cellRegionMaster.size(), label(-1));

    // Get point-cell addressing
    const labelListList& pc = mesh_.pointCells();

    forAll (splitPointsToUnrefine, i)
    {
        const labelList& curPc = pc[splitPointsToUnrefine[i]];

        label curRegion = -1;

        forAll (curPc, curPcI)
        {
            if (curRegion == -1)
            {
                // First region found.  Grab it
                curRegion = cellRegion[curPc[curPcI]];
            }
            else
            {
                // Region already found.  Check that all other cells that
                // touch this point have the same region
                if (curRegion != cellRegion[curPc[curPcI]])
                {
                    // Error: different region cells touching in split point
                    // This is not a valid unrefinement pattern
                    FatalErrorInFunction
                        << "Different region cells touching in split point."
                        << abort(FatalError);
                }
            }
        }

        // Record point region master
        if (curRegion > -1)
        {
            pointRegionMaster[curRegion] = splitPointsToUnrefine[i];
        }
        else
        {
            // Error: Cannot find region for point
            FatalErrorInFunction
                << "Different region cells touching in split point."
                << abort(FatalError);
        }
    }

    // Insert all commands to combine cells
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
//         pointRegionMaster,
        cellRegionMaster,
        meshMod
    );

    // Remove the 8 cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitPointsToUnrefine, i)
    {
        label pointi = splitPointsToUnrefine[i];

        const labelList& pCells = mesh_.pointCells(pointi);

        label masterCelli = min(pCells);

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyhedralRefinement::polyhedralRefinement
(
    const polyMesh& mesh,
    const dictionary& dict,
    const bool read
)
:
    refinement(mesh, dict, read)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyhedralRefinement::~polyhedralRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::polyhedralRefinement::consistentUnrefinement
(
    const labelList& unrefinementPointCandidates,
    const bool maxSet
) const
{
    if (debug)
    {
        InfoInFunction<< "Setting split points to unrefine." << endl;
    }

    // Get necessary mesh data
    const label nPoints = mesh_.nPoints();

    // PART 1: Mark all split points in the mesh (points that can be unrefined)
    boolList splitPointsMarkup(nPoints, false);

    // Algorithm: split point is uniquely defined as a point that:
    // 1. Has pointLevel_ > 0 (obviously),
    // 2. A point that has the same pointLevel_ as ALL of the points of its
    //    faces. In other words, for each point, we will look through all the
    //    faces of the point. For each face, we will visit points and check the
    //    point level of all of these points. All point levels must be the same
    //    for this point candidate to be split point. This is quite useful since
    //    there is no need to store the refinement history

    // Get necessary mesh data
    const faceList& meshFaces = mesh_.faces();
    const labelListList& meshPointFaces = mesh_.pointFaces();

    // Loop through all points
    forAll (meshPointFaces, pointI)
    {
        // Get point level of this point
        const label& centralPointLevel = pointLevel_[pointI];

        if (centralPointLevel < 1)
        {
            // Point can't be unrefined as its level is either 0 or
            // invalid. Continue immediately
            continue;
        }

        // Flag to see whether this is a split point candidate
        bool splitPointCandidate = true;

        // Get face labels for this point
        const labelList& pFaces = meshPointFaces[pointI];

        // Loop through all point faces
        forAll (pFaces, i)
        {
            // Get face index and the face
            const label& faceI = pFaces[i];
            const face& curFace = meshFaces[faceI];

            // Loop through points of the face
            forAll (curFace, j)
            {
                // Get point index
                const label& pointJ = curFace[j];

                if (pointLevel_[pointJ] != centralPointLevel)
                {
                    // Point levels are different, this can't be a split point,
                    // set flag to false and break immediatelly
                    splitPointCandidate = false;
                    break;
                }
                // else: this is still potential split point candidate so
                //       there's nothing to do
            } // End for all points of this face

            // Check whether this can't be a split point already and break out
            // immediately
            if (!splitPointCandidate)
            {
                break;
            }
        } // End for all point faces

        // At this point, if the flag is still true, this is a split point
        if (splitPointCandidate)
        {
            splitPointsMarkup[pointI] = true;
        }
    }

    // Note: if there is no dynamic load balancing, points at the boundary
    // cannot be split points by definition. However, in dynamic load balancing
    // runs, it is possible that a split point ends on processor boundary, in
    // which case we will simply avoid (actually delay) unrefining until this
    // becomes internal point again. VV, 4/Jul/2018.
    const label nInternalFaces = mesh_.nInternalFaces();
    const label nFaces = mesh_.nFaces();

    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
    {
        // Get the face and make sure that the points are unmarked
        const face& f = meshFaces[faceI];

        forAll (f, fpI)
        {
            splitPointsMarkup[f[fpI]] = false;
        }
    }

    // PART 2: Mark all unrefinement point candidates that are split points at
    // the same time (basically the intersection of split points and candidates)

    // Create markup field of split points to unrefine
    // True: this is a split point which should be unrefined
    // False: this is either not a split point or it shouldn't be unrefined
    PackedBoolList splitPointsToUnrefine(nPoints);

    // Loop through all unrefinement candidates
    forAll (unrefinementPointCandidates, i)
    {
        // Get point index
        const label& pointI = unrefinementPointCandidates[i];

        if (splitPointsMarkup[pointI])
        {
            // This is a split point, mark it for unrefinement
            splitPointsToUnrefine.set(pointI);
        }
    }


    // PART 3: Ensure face consistent (2:1 constraint) and possibly point
    // consistent (4:1 constraint) unrefinement

    // Get necessary mesh data
    const label nCells = mesh_.nCells();
    const labelListList& meshPointCells = mesh_.pointCells();

    // Count number of removed cells from unrefinement (cells that will not be
    // unrefined) in each iteration and number of iterations
    label nChanged = 0;

    while (true)
    {
        // First, create cells to unrefine (all cells sharing point to unrefine)
        PackedBoolList cellsToUnrefine(nCells, false);

        // Loop through all split points to unrefine
        forAll (splitPointsToUnrefine, pointI)
        {
            if (splitPointsToUnrefine.get(pointI))
            {
                // This split point is marked for unrefinement, collect all of
                // its cells
                const labelList& pCells = meshPointCells[pointI];
                forAll (pCells, i)
                {
                    cellsToUnrefine.set(pCells[i]);
                }
            }
        }

        // Reset number of removed cells from unrefinement for this iteration
        nChanged = faceConsistentUnrefinement(cellsToUnrefine, maxSet);

        if (edgeBasedConsistency_)
        {
            // Check for 4:1 edge based consistent unrefinement. Updates
            // cellsToUnrefine and returns number of removed cells from
            // unrefinement in this iteration
            nChanged += edgeConsistentUnrefinement(cellsToUnrefine, maxSet);
        }

        // Global reduction
        reduce(nChanged, sumOp<label>());

        // If we have removed at least one cell from unrefinement, we need to
        // protect its split points as well from unrefinement
        if (nChanged == 0)
        {
            break;
        }

        // Get point cells
        const labelListList& meshPointCells = mesh_.pointCells();

        // Loop through all split points to unrefine
        forAll (splitPointsToUnrefine, pointI)
        {
            if (splitPointsToUnrefine.get(pointI))
            {
                // This is a split point for unrefinement, get the cells
                const labelList& pCells = meshPointCells[pointI];

                // Loop through all point cells
                forAll (pCells, i)
                {
                    if (!cellsToUnrefine[pCells[i]])
                    {
                        // Cell must not be refined, remove point from
                        // unrefinement as well
                        splitPointsToUnrefine.unset(pointI);
                        break;
                    }
                }
            }
        }

        if (debug)
        {
            Pout<< "Removed " << nChanged << " cells" << endl;
        }

    }

    // Convert back to labelList.
    label nSet = splitPointsToUnrefine.count();

    labelList newPointsToUnrefine(nSet);
    nSet = 0;

    forAll(splitPointsToUnrefine, pointI)
    {
        if (splitPointsToUnrefine.get(pointI))
        {
            newPointsToUnrefine[nSet++] = pointI;
        }
    }

    return newPointsToUnrefine;
}


// ************************************************************************* //
