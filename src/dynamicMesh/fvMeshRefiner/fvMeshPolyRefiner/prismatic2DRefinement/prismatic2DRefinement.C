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

#include "prismatic2DRefinement.H"
#include "polyTopoChange.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "syncTools.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"
#include "foamMeshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(prismatic2DRefinement, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::prismatic2DRefinement::getAnchorLevel
(
    const label faceI,
    const label nPoints
) const
{
    // Get the face
    const face& f = mesh_.faces()[faceI];

    // Sanity check for expected number of points
    if (nPoints != 3 && nPoints != 4)
    {
        FatalErrorInFunction
            << "Trying to find anchor level with " << nPoints << " points"
            << " smaller than anchor level." << nl
            << "Only nPoints = 3 and 4 are supported." << endl
            << abort(FatalError);
    }

    // Sanity check: if we are expecting to find 4 points on a face which is not
    // on special patch (empty or wedge) and we find a face with 3 points, issue
    // an error
    if (f.size() <= 3 && nPoints == 4)
    {
        FatalErrorInFunction
            << "Expected to find at least 4 points with level lower than"
            << " anchor level." << nl
            << "Make sure to call this function with nPoints = 4 only if"
            << " you do not expect triangular faces." << endl
            << abort(FatalError);
    }

    if (f.size() <= nPoints)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        const label& ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) >= nPoints)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel + 1) >= nPoints)
        {
            return ownLevel + 1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::prismatic2DRefinement::appendFaceSplitInfo
(
    const label& faceI,
    const boolList& edgeOnPatchToCut,
    const labelList& edgeMidPoint,
    DynamicList<label>& splitFacesIntoTwo,
    DynamicList<Pair<label> >& splitFacesEmptyEdges
) const
{
    // First append the face into the list
    splitFacesIntoTwo.append(faceI);

    // Grab all edges of the face
    const labelList& curEdges = mesh_.faceEdges()[faceI];

    // Create placeholders for two edge levels. Initialise with -1
    // for sanity checks later on
    label specialPatchEdgeI = -1;
    label specialPatchEdgeJ = -1;

    // Count number of edges for sanity checks
    label nEdgesOnSpecialPatch = 0;

    // Collect the two edge labels found on the special patch (empty or wedge)
    forAll (curEdges, fi)
    {
        // Get edge index
        const label edgeI = curEdges[fi];

        if (edgeOnPatchToCut[edgeI])
        {
            // This edge is on special patch (empty or wedge), check whether it
            // is first, second or invalid
            switch (nEdgesOnSpecialPatch++)
            {
                case 0:
                    specialPatchEdgeI = edgeI;
                    break;
                case 1:
                    specialPatchEdgeJ = edgeI;
                    break;
                default:
                    FatalErrorInFunction
                        << "Found more than two edges on face " << faceI
                        << " on the special patch (empty or wedge)." << nl
                        << "Either this is not a valid 2D mesh or"
                        << " we are visiting wrong faces." << endl
                        << abort(FatalError);
            }

        } // End if edge on special patch (empty or wedge)

    } // End loop over all edges

    // Debug: additional check whether the two edges are marked for
    // refinement (they should be)
    if (debug)
    {
        if (edgeMidPoint[specialPatchEdgeI] == -1)
        {
            FatalErrorInFunction
                << "Empty patch edge with index: " << specialPatchEdgeI
                << " not marked for splitting" << nl
                << "Check edgeMidPoint selection algorithm." << endl
                << abort(FatalError);
        }

        if (edgeMidPoint[specialPatchEdgeJ] == -1)
        {
            FatalErrorInFunction
                << "Empty patch edge with index: " << specialPatchEdgeJ
                << " not marked for splitting" << nl
                << "Check edgeMidPoint selection algorithm." << endl
                << abort(FatalError);
        }
    }

    // At this point, we should have the two edges we were looking
    // for, collect them into the list with additional sanity check
    if (specialPatchEdgeI > -1 && specialPatchEdgeJ > -1)
    {
        // Append the list of two edges in increasing order (just in
        // case we end up needing this information
        if (specialPatchEdgeI < specialPatchEdgeJ)
        {
            splitFacesEmptyEdges.append
            (
                Pair<label>(specialPatchEdgeI, specialPatchEdgeJ)
            );
        }
        else
        {
            splitFacesEmptyEdges.append
            (
                Pair<label>(specialPatchEdgeJ, specialPatchEdgeI)
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Found invalid indices for edges on special patches:" << nl
            << "specialPatchEdgeI: " << specialPatchEdgeI
            << ", specialPatchEdgeJ: " << specialPatchEdgeJ << nl
            << "Something went wrong. Check face edges." << endl
            << abort(FatalError);
    }
}


void Foam::prismatic2DRefinement::setNewFaceNeighbours
(
    const HashTable
    <
        label,
        Pair<label>,
        Hash<FixedList<label, 2> >
    >& pointCellToAddedCellMap,
    const labelListList& cellAddedCells,
    const label& faceI,
    const label& pointI,

    label& own,
    label& nei
) const
{
    // Get anchor cell for this anchor point on owner side
    const label ownCellI = mesh_.faceOwner()[faceI];

    // Get cell added cells for owner cell
    const labelList& cAddedOwn = cellAddedCells[ownCellI];

    // If the cell added cells list is not empty, return the necessary cell
    // index fetched with pointCellToAddedCellMap
    const Pair<label> pointOwnerCellPair(pointI, ownCellI);

    if (cAddedOwn.empty())
    {
        // Cell added cells is empty for owner, fetch original owner
        own = ownCellI;
    }
    else if (!pointCellToAddedCellMap.found(pointOwnerCellPair))
    {
        // Point-cell pair not found, meaning that we need to search for the
        // closest point which has lower level than this point. This happens if
        // we are splitting a face that has already been split but has different
        // owner/neighbour levels (e.g. owner was refined, but the neighbour was
        // not in the previous time step). It is possible that we end up with
        // splitting this face with point levels e.g. (1 1 0 0). Therefore, we
        // need to search for the point with minimum level on the edges sharing
        // this point.

        // Find point in the local faces addressing
        const face& f = mesh_.faces()[faceI];
        const label fpI = findIndex(f, pointI);

        if (fpI == -1)
        {
            FatalErrorInFunction
                << "Point: " << pointI << " not found in face: " << f
                << ", with face index: " << faceI << nl
                << "Point must belong to the face." << nl
                << "Error encounted when dealing with owner cell: "
                << ownCellI << endl
                << abort(FatalError);
        }

        // Find the global point index with minimum edge connected level
        const label anchorPointI = findMinEdgeConnectedLevel
        (
            fpI,                      // current face point
            faceI,                    // face index
            f,                        // face
            mesh_.faceEdges()[faceI], // face edges
            mesh_.edges()             // mesh edges
        );

        // If the point is the same as pointI, we did not find any valid point
        if (anchorPointI == pointI)
        {
            FatalErrorInFunction
                << "Could not find different adjacent anchor point." << nl
                << "pointI: " << pointI << " faceI: " << faceI << nl
                << "Error encounted when dealing with owner cell: "
                << ownCellI << endl
                << abort(FatalError);
        }

        // Set owner cell
        own = cAddedOwn
        [
            pointCellToAddedCellMap[Pair<label>(anchorPointI, ownCellI)]
        ];
    }
    else
    {
        // Cell added cells is not empty and the mapping is found. Set owner
        // from mapping
        own = cAddedOwn[pointCellToAddedCellMap[pointOwnerCellPair]];
    }

    if (mesh_.isInternalFace(faceI))
    {
        // Get anchor cell for this anchor point on neighbour side
        const label neiCellI = mesh_.faceNeighbour()[faceI];

        // Get cell added cells for neighbour cell
        const labelList& cAddedNei = cellAddedCells[neiCellI];

        // If the cell added cells list is not empty, return the necessary cell
        // index fetched with pointCellToAddedCellMap

        const Pair<label> pointNeighbourCellPair(pointI, neiCellI);

        if (cAddedNei.empty())
        {
            // Cell added cells is empty for neighbour, fetch original neighbour
            nei = neiCellI;
        }
        else if (!pointCellToAddedCellMap.found(pointNeighbourCellPair))
        {
            // Point-cell pair not found, meaning that we need to search for the
            // closest point which has lower level than this point. This happens
            // if we are splitting a face that has already been split but has
            // different owner/neighbour levels (e.g. owner was refined, but the
            // neighbour was not in the previous time step). It is possible that
            // we end up with splitting this face with point levels e.g. (1 1 0
            // 0). Therefore, we need to search for the point with minimum level
            // on the edges sharing this point.

            // Find point in the local faces addressing
            const face& f = mesh_.faces()[faceI];
            const label fpI = findIndex(f, pointI);

            if (fpI == -1)
            {
                FatalErrorInFunction
                    << "Point: " << pointI << " not found in face: " << f
                    << ", with face index: " << faceI << nl
                    << "Point must belong to the face." << nl
                    << "Error encounted when dealing with neighbour cell: "
                    << neiCellI << endl
                    << abort(FatalError);
            }

            // Find the global point index with minimum edge connected level
            const label anchorPointI = findMinEdgeConnectedLevel
            (
                fpI,                      // current face point
                faceI,                    // face index
                f,                        // face
                mesh_.faceEdges()[faceI], // face edges
                mesh_.edges()             // mesh edges
            );

            // If the point is the same as pointI, we did not find any valid point
            if (anchorPointI == pointI)
            {
                FatalErrorInFunction
                    << "Could not find different adjacent anchor point."<< nl
                    << "pointI: " << pointI << " faceI: " << faceI
                    << "Error encounted when dealing with neighbour cell: "
                    << neiCellI << endl
                    << abort(FatalError);
            }

            // Set owner cell
            nei = cAddedNei
            [
                pointCellToAddedCellMap[Pair<label>(anchorPointI, neiCellI)]
            ];
        }
        else
        {
            // Cell added cells is not empty and the mapping is found. Set
            // neighbour from mapping
            nei = cAddedNei[pointCellToAddedCellMap[pointNeighbourCellPair]];
        }
    }
    else
    {
        // Boundary face: set neighbour to -1
        nei = -1;
    }
}


Foam::label Foam::prismatic2DRefinement::findMinEdgeConnectedLevel
(
    const label& fpI,
    const label& faceI,
    const face& f,
    const labelList& fEdges,
    const edgeList& meshEdges
) const
{
    // Get point index and initialize anchor point
    const label& pointI = f[fpI];

    label anchorPointI = pointI;

    // Check the other point on edge starting with pointI
    const label& edgeIndexAfterPoint = fEdges[fpI];
    const edge& edgeAfter = meshEdges[edgeIndexAfterPoint];

    // Get other point on edge after
    const label& pointAfter = edgeAfter.otherVertex(pointI);

    if (pointLevel_[pointAfter] < pointLevel_[anchorPointI])
    {
        anchorPointI = pointAfter;
    }

    // Check the other point on edge ending with pointI
    const label& edgeIndexBeforePoint = fEdges[f.rcIndex(fpI)];
    const edge& edgeBefore = meshEdges[edgeIndexBeforePoint];

    // Get other point on edge before
    const label& pointBefore = edgeBefore.otherVertex(pointI);

    if (pointLevel_[pointBefore] < pointLevel_[anchorPointI])
    {
        anchorPointI = pointBefore;
    }

    return anchorPointI;
}


void Foam::prismatic2DRefinement::addFaceMids
(
    const labelList& faceMidPoint,
    const boolList& faceOnPatchToCut,
    const label& faceI,
    const label& cellI,
    face& newFace
) const
{
    // b) faceMidPoint for this face
    newFace[1] = faceMidPoint[faceI];

    // c) faceMidPoint for the face on the other side
    // The other face is uniquely defined as the other face of the same cell
    // which is on special patch (empty or wedge)

    // Get the cell
    const cell& cFaces = mesh_.cells()[cellI];

    // Loop through cell faces
    forAll(cFaces, i)
    {
        // Get face index
        const label& faceJ = cFaces[i];

        if (faceOnPatchToCut[faceJ] && (faceI != faceJ))
        {
            // This is the face we're looking for, add its
            // midpoint and double check if it is valid
            if (faceMidPoint[faceJ] > -1)
            {
                newFace[2] = faceMidPoint[faceJ];
                break;
            }
            else
            {
                FatalErrorInFunction
                    << "Other face: " << faceJ
                    << " has not been selected for splitting,"
                    << " while the face on original side: "<< faceI
                    <<" has been selected." << endl
                    << abort(FatalError);
            }
        } // End if this is our face
    } // End for all cell faces
}


void Foam::prismatic2DRefinement::checkNewFaceOrientation
(
    polyTopoChange& meshMod,
    const label& faceI,
    const face& newFace
) const
{
    // Get mesh cell centres
    const vectorField& meshCellCentres = mesh_.cellCentres();

    if (mesh_.isInternalFace(faceI))
    {
        // Get old owner/neighbour indices
        const label oldOwn = mesh_.faceOwner()[faceI];
        const label oldNei = mesh_.faceNeighbour()[faceI];

        // Print info only with deep debug level
        if (debug > 1)
        {
            Pout<< "Split infternal face: " << faceI
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
        // Get face centres and old owner
        const vectorField& meshFaceCentres = mesh_.faceCentres();
        const label oldOwn = mesh_.faceOwner()[faceI];

        // Print info only with deep debug level
        if (debug > 1)
        {
            Pout<< "Split boundary face: " << faceI
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
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::prismatic2DRefinement::setRefinement
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

    // PART 1: Mark cells for refinement

    // Bool list that marks cells which will be refined
    PackedBoolList refineCellsMask(mesh_.nCells(), false);
    forAll(cellsToRefine, i)
    {
        // Simply mark the cell as refined, there are no additional points to
        // add (in cell centre for example)
        refineCellsMask.set(cellsToRefine[i]);
    }

    // Write out cells to refine as a cell set for debug
    if (debug)
    {
        // Note: cellSet is actually a hash table of labels
        cellSet splitCells(mesh_, "splitCells", cellsToRefine.size());

        forAll(refineCellsMask, cellI)
        {
            if (refineCellsMask.get(cellI))
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

    // First mark faces and edges on special patches (empty or wedge). This data
    // is used in PART 2 (collecting edges) and also PART 3 (collecting faces)
    boolList faceOnPatchToCut(mesh_.nFaces(), false);
    boolList edgeOnPatchToCut(mesh_.nEdges(), false);

    // Get boundary
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();

    // Get face-edge addressing (for each face, a list of edges)
    const labelListList& meshFaceEdges = mesh_.faceEdges();

    // Loop through all patches
    forAll(boundaryMesh, patchI)
    {
        // Get current patch
        const polyPatch& curPatch = boundaryMesh[patchI];

        // Check whether this patch is special (empty or wedge)
        if (isA<emptyPolyPatch>(curPatch) || isA<wedgePolyPatch>(curPatch))
        {
            // Get start and end face labels
            const label startFaceI = curPatch.start();
            const label endFaceI = startFaceI + curPatch.size();

            // Mark all the faces and edges on the patch
            for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
            {
                // Mark face
                faceOnPatchToCut[faceI] = true;

                // Get edges of this face
                const labelList& curEdges = meshFaceEdges[faceI];

                // Mark all edges
                forAll(curEdges, i)
                {
                    edgeOnPatchToCut[curEdges[i]] = true;
                }
            }
        }
    }

    // Now that we have marked faces and edges on special patches (empty or
    // wedge), let's collect refined edges. Refined edges are defined by having
    // both their point levels <= cell level, i.e. if any cell that gets split
    // uses this edge and the edge is on special patch (empty or wedge), the
    // edge needs to be split

    // Get necessary mesh data
    const labelListList& meshCellEdges = mesh_.cellEdges();
    const edgeList& meshEdges = mesh_.edges();

    // Mid points for refined edge:
    // No need to split edge = -1
    // Label of introduced mid point > -1
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over refined cells
    forAll(cellsToRefine, i)
    {
        // Get cell index
        const label& cellI = cellsToRefine[i];

        // Get edges of this cell
        const labelList& cEdges = meshCellEdges[cellI];

        forAll (cEdges, j)
        {
            // Get edge index and edge
            const label& edgeI = cEdges[j];
            const edge& e = meshEdges[edgeI];

            if
            (
                edgeOnPatchToCut[edgeI]
             && pointLevel_[e[0]] <= cellLevel_[cellI]
             && pointLevel_[e[1]] <= cellLevel_[cellI]
            )
            {
                // Point levels of both edge points are <= cell level, mark
                // edge for splitting
                edgeMidPoint[edgeI] = 12345;
            }
            // Else nothing to do: can't split edges that are not on special
            // patch (empty or wedge)
        } // End for all edges of the refined cell
    } // End for all cells

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

    // Now that the refinement trigger is synced, introduce edge points

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
                        e[0],              // Appended point, no master ID
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


    // PART 3: Collect faces for refinement. Faces need to be collected in two
    // distinct categories:
    // 1. Faces found on special patches (empty or wedge) that will be split
    //    into n faces (where n is the number of edges per face),
    // 2. Faces not on special patches (empty or wedge) that will be always
    //    split into two faces. For each of these faces, collect the two edges
    //    found on opposing sides of the special patch (empty or wedge).

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Allocating face midpoints and collecting faces that are"
            << " not on special patch (empty or wedge)."
            << endl;
    }

    // Get face anchor level based on the face type. For split face found on
    // special patch (empty or wedge), it is guaranteeed that we will have at
    // least 3 points with level <= anchor level. For split face not on special
    // (empty or wedge patch), it is guaranteed that we will have at least 4
    // points with level <= anchor level. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());
    for (label faceI = 0; faceI < mesh_.nFaces(); ++faceI)
    {
        if (faceOnPatchToCut[faceI])
        {
            // Face on special patch, at least 3 points need to have
            // level <= anchor level
            faceAnchorLevel[faceI] = getAnchorLevel(faceI, 3);
        }
        else
        {
            // Face not on special patch, at least 4 points need to have
            // level <= anchor level
            faceAnchorLevel[faceI] = getAnchorLevel(faceI, 4);
        }
    }

    // Split faces on special patches (empty or wedge) will be collected in
    // faceMidPoint list:
    // Not refined = -1
    // Shall be refined > -1 (label of added mid point)
    labelList faceMidPoint(mesh_.nFaces(), -1);

    // Split faces not on special patches (empty or wedge) will be collected
    // into splitFacesIntoTwo dynamic list. For each of these faces, we also
    // need to collect its two edges that are found on special patch (empty or
    // wedge)

    // Allocate enough storage to prevent excessive resizing
    const label nSplitFacesIntoTwo = 3*cellsToRefine.size();
    DynamicList<label> splitFacesIntoTwo(nSplitFacesIntoTwo);
    DynamicList<Pair<label> > splitFacesEmptyEdges(nSplitFacesIntoTwo);

    // Get necessary mesh data
    const labelList& meshFaceOwner = mesh_.faceOwner();
    const labelList& meshFaceNeighbour = mesh_.faceNeighbour();

    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Internal faces: look at cells on both sides. Uniquely determined since
    // the face itself is guaranteed to be same level as most refined neighbour
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        // Get owner data
        const label& own = meshFaceOwner[faceI];
        const label& ownLevel = cellLevel_[own];
        const label newOwnLevel = ownLevel + refineCellsMask.get(own);

        // Get neighbour data
        const label& nei = meshFaceNeighbour[faceI];
        const label& neiLevel = cellLevel_[nei];
        const label newNeiLevel = neiLevel + refineCellsMask.get(nei);

        if
        (
            newOwnLevel > faceAnchorLevel[faceI]
         || newNeiLevel > faceAnchorLevel[faceI]
        )
        {
            // Note: this is internal face so we don't need to check whether the
            // face is on special patch (empty or wedge). It can't be by
            // definition

            // Does two things:
            // 1. Appends the face to splitFacesIntoTwo list
            // 2. Append the two edges on special patch (empty or wedge) to
            //    splitFaceEmptyEdges list
            appendFaceSplitInfo
            (
                faceI,
                edgeOnPatchToCut,
                edgeMidPoint,
                splitFacesIntoTwo,
                splitFacesEmptyEdges
            );
        } // End whether the face needs to be considered (split)
    } // End loop over all internal faces

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
            const label newOwnLevel = ownLevel + refineCellsMask.get(own);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap the list which now contains data from the other side
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        forAll(newNeiLevel, i)
        {
            // Get face index
            const label faceI = i + nInternalFaces;

            // Get owner data (neighbour is available from before)
            const label& own = meshFaceOwner[faceI];
            const label& ownLevel = cellLevel_[own];
            const label newOwnLevel = ownLevel + refineCellsMask.get(own);

            if
            (
                newOwnLevel > faceAnchorLevel[faceI]
             || newNeiLevel[i] > faceAnchorLevel[faceI]
            )
            {
                if (faceOnPatchToCut[faceI])
                {
                    // This face is on the special patch (empty or wedge) and
                    // will be split into n faces (n is the number of edges for
                    // this face) and the face mid point will be added. Mark the
                    // face for splitting
                    faceMidPoint[faceI] = 12345;
                }
                else
                {
                    // Does two things:
                    // 1. Appends the face to splitFacesIntoTwo list
                    // 2. Append the two edges on special patch (empty or wedge)
                    //    to splitFaceEmptyEdges list
                    appendFaceSplitInfo
                    (
                        faceI,
                        edgeOnPatchToCut,
                        edgeMidPoint,
                        splitFacesIntoTwo,
                        splitFacesEmptyEdges
                    );
                } // End if the face is not on special patch (empty or wedge)
            } // End whether the face needs to be considered
        } // End loop over all boundary faces
    } // End memory management for syncing owner/neighbour face levels


    // Add face points. Note: no need to sync face mid points (as we did for
    // edge mid points) since processor faces do not introduce new points, only
    // faces on special patch (empty or wedge) do

    // Get face centres
    const vectorField& meshFaceCentres = mesh_.faceCentres();

    // Loop through faces on special patches (empty or wedge) only
    forAll (boundaryMesh, patchI)
    {
        // Get current patch
        const polyPatch& curPatch = boundaryMesh[patchI];

        // Check whether this patch is special (empty or wedge)
        if (isA<emptyPolyPatch>(curPatch) || isA<wedgePolyPatch>(curPatch))
        {
            // Get start and face labels
            const label startFaceI = curPatch.start();
            const label endFaceI = startFaceI + curPatch.size();

            // Loop through all special patch faces (global indexing)
            for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
            {
                if (faceMidPoint[faceI] > -1)
                {
                    const face& f = mesh_.faces()[faceI];

                    // Face on special patch (empty or wedge) marked to be
                    // split. Add the point at face centre and replace
                    // faceMidPoint with new point label

                    faceMidPoint[faceI] = meshMod.setAction
                    (
                        polyAddPoint
                        (
                            meshFaceCentres[faceI], // Point
                            f[0],                   // No master ID
                            -1,                     // Zone for point
                            true                    // Supports a cell
                        )
                    );

                    // Determine the level of the corner points and midpoint will
                    // be one higher.
                    newPointLevel(faceMidPoint[faceI]) = faceAnchorLevel[faceI]+1;
                } // End if face marked for splitting
            } // End loop over all faces on special patch (empty or wedge)
        } // End if special patch check
    } // End loop for all patches

    // Write out all split faces as a face set for debugging
    if (debug)
    {
        // Create faceSet containing all faces that need to be split into n
        // faces (n is the number of edges on the face)
        faceSet splitNFaces(mesh_, "splitNFaces", 3*cellsToRefine.size());

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                splitNFaces.insert(faceI);
            }
        }

        Pout<< FUNCTION_NAME << nl
            << "Writing " << splitNFaces.size()
            << " faces to split in N to faceSet " << splitNFaces.objectPath()
            << endl;

        splitNFaces.write();

        // Create faceSet containing all faces that need to be split into 2
        // faces (faces not on special patch: empty or wedge)
        faceSet splitTwoFaces(mesh_, "splitTwoFaces", 3*cellsToRefine.size());

        forAll (splitFacesIntoTwo, i)
        {
            // Insert face index into splitTwoFaces
            splitTwoFaces.insert(splitFacesIntoTwo[i]);
        }

        Pout<< FUNCTION_NAME << nl
            << "Writing " << splitTwoFaces.size()
            << " faces to split in two to faceSet "
            << splitTwoFaces.objectPath() << endl;

        splitTwoFaces.write();
    }


    // Now we have all the information we need to perform the refinement and we
    // no longer need to refer to cellsToRefine_. The information is in:
    // - refineCellsMask = true : cell needs to be split
    // - edgeMidPoint >= 0     : edge on special patch that needs to be split
    // - faceMidPoint >= 0     : face on special patch that needs to be split
    //                           into n faces (where n is the number of edges)
    // - splitFacesIntoTwo     : list of faces that need to be split into two
    //                           (face not on special patch)
    // - splitFacesEmptyEdges  : holds the two edges of the face which needs to
    //                           be split into two


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
    forAll(refineCellsMask, cellI)
    {
        if (refineCellsMask.get(cellI))
        {
            // The cell will be refined, set capacity to 8 to prevent excessive
            // resizing
            cellAnchorPointsDynamic[cellI].setCapacity(8);
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
                refineCellsMask.get(cellI)
             && pointLevel_[pointI] <= cellLevel_[cellI]
            )
            {
                // This point cell is marked for refinement and its point level
                // is smaller or equal to cell level, append the point
                cellAnchorPointsDynamic[cellI].append(pointI);
            }
        }
    }

    // Loop through all cells and check whether at least 6 anchor points
    // have been found (minimum requirement for a triangular prism)

    // Collect cellAnchorPoint into a List<labelList> instead of
    // List<dynamicList>
    labelListList cellAnchorPoints(mesh_.nCells());

    // Get cell points for error output
    const labelListList& meshCellPoints = mesh_.cellPoints();

    forAll(refineCellsMask, cellI)
    {
        // First some sanity checks
        if (refineCellsMask.get(cellI))
        {
            // Cell selected for refinement
            if (cellAnchorPointsDynamic[cellI].size() < 6)
            {
                // Cell has less than 6 anchor points. Issue an error and report
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
            else if (cellAnchorPointsDynamic[cellI].size() % 2 != 0)
            {
                // Cell has odd number of anchor points. This is not allowed and
                // indicates an invalid mesh
                const labelList& cPoints = meshCellPoints[cellI];

                FatalErrorInFunction
                    << "Cell " << cellI
                    << " of level " << cellLevel_[cellI]
                    << " has odd number of anchor points"
                    << " (should be even for 2D mesh). "
                    << "cellPoints:" << cPoints << endl
                    << "pointLevels:"
                    << IndirectList<label>(pointLevel_, cPoints)() << endl
                    << abort(FatalError);
            }
        }

        // Tranfer the dynamic list for each cell into an ordinary list
        cellAnchorPoints[cellI].transfer(cellAnchorPointsDynamic[cellI]);
    }


    // PART 5: Add the cells

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << " Adding cells."
            << endl;
    }

    // We should have exactly n new cells per each split cell, where n is the
    // number of anchor points in a cell divided by two. In order to determine
    // owner/neighbours of new and modified faces, we need to know which cell
    // came from which point. The mapping is not uniquely defined as in
    // polyhedralRefinement when we had 1 point = 1 cell. Here, we have two
    // points that correspond to a single cell, one on one side of the special
    // patch and the other on other side. This information will be collected in
    // a HashTable<label, Pair<label> >, where the key will be a pair of
    // global point index and global cell index, while the value is local index
    // into cellAddedCells list
    labelListList cellAddedCells(mesh_.nCells());
    HashTable<label, Pair<label>, Hash<FixedList<label, 2> > >
        pointCellToAddedCellMap(6*cellsToRefine.size());

    // Get mesh data
    const meshCellZones& cellZones = mesh_.cellZones();
    const cellList& meshCells = mesh_.cells();
    const faceList& meshFaces = mesh_.faces();
    const labelListList& meshPointEdges = mesh_.pointEdges();


    // Loop through all faces
    forAll(faceMidPoint, faceI)
    {
        // Check whether this face needs to be split into n faces
        if (faceMidPoint[faceI] > -1)
        {
            // This is a face that will be split and is on special patch (empty
            // or wedge) by definition, get the cell index by looking at owner
            // only
            const label& cellI = meshFaceOwner[faceI];

            // Get cell added cells
            labelList& cAdded = cellAddedCells[cellI];

            // Check whether the cell added cells are empty. This means that we
            // haven't visited the first face yet. If it is not empty, we have
            // already visited one face, which is enough
            if (cAdded.empty())
            {
                // First face that hasn't been visited. Start adding cells
                // point-by-point and keep track of mapping necessary for
                // splitting other (not on special patch: empty or wedge) faces
                // into two

                // Set the total number of added cells to number of anchors
                // divided by two. Note: number of anchors needs to be an even
                // number (6 for triangular prism, 8 for hex, 10 for pentagonal
                // prism, etc.)
                const labelList& cAnchors = cellAnchorPoints[cellI];
                cAdded.setSize(cAnchors.size()/2);

                // Helper variable to distinguish between first and successive
                // cells (first will have the original index)
                label cellCounter = 0;

                // Get current face
                const face& f = meshFaces[faceI];

                // Loop through face points
                forAll (f, fpI)
                {
                    // Get point index
                    const label& pointI = f[fpI];

                    // Find anchor point in local list if present
                    const label anchorI = findIndex(cAnchors, pointI);

                    if (anchorI != -1)
                    {
                        // This point is anchor, add the cell

                        if (cellCounter == 0)
                        {
                            // This is first cell, simply set the existing index
                            cAdded[cellCounter] = cellI;
                        }
                        else
                        {
                            // Other cells, need to add the cells
                            cAdded[cellCounter] = meshMod.setAction
                            (
                                polyAddCell
                                (
                                    -1,                         // M. point
                                    -1,                         // M. edge
                                    -1,                         // M. face
                                    cellI,                      // M. cell
                                    cellZones.whichZone(cellI)  // M. zone
                                )
                            );
                        }

                        // Update cell level
                        newCellLevel(cAdded[cellCounter]) = cellLevel_[cellI]
                        + 1;

                        // Collect the point-cell mapping into local index
                        // of cell added cells for point on this side
                        pointCellToAddedCellMap.insert
                        (
                            Pair<label>(pointI, cellI),
                            cellCounter
                        );

                        // This is only one side, we need to also collect
                        // the other point on the other side. Get point edges
                        // for this point
                        const labelList& pEdges =
                            meshPointEdges[pointI];

                        // Loop through point edges
                        forAll (pEdges, peI)
                        {
                            // Get the edge index
                            const label& edgeI = pEdges[peI];

                            if (!edgeOnPatchToCut[edgeI])
                            {
                                // Edge is not on special patch (empty or
                                // wedge), therefore this is the edge we're
                                // looking for. Collect the other point of the
                                // edge
                                const label pointJ =
                                    meshEdges[edgeI].otherVertex(pointI);

                                // Collect the point-cell mapping into local
                                // index of cell added cells for this point
                                pointCellToAddedCellMap.insert
                                (
                                    Pair<label>(pointJ, cellI),
                                    cellCounter
                                );
                            }
                        }

                        // Now we have finished adding the cell and also
                        // adding the necessary mapping for this added cell

                        // Increment cell counter
                        ++cellCounter;

                    } // Else point is not anchor: nothing to do
                } // End loop over all face points

                // Sanity check: number of counted cells must be
                // equal to size of cellAddedCells. This means that
                // we have correctly marked the anchor points
                if (cellCounter != cAdded.size())
                {
                    FatalErrorInFunction
                        << "Problem while adding cells."
                        << nl
                        << "Going through base face on special patch"
                        << " (empty or wedge) and adding cells, we collected: "
                        << cellCounter << " cells."
                        << nl
                        << "While the number of anchor points is: "
                        << cAnchors.size()
                        << nl
                        << "The number of added cells based on number of anchor"
                        << " points is: "
                        << cAdded.size()
                        << nl
                        << "Additional information: "
                        << nl
                        << "cellI: " << cellI << ", faceI: " << faceI
                        << abort(FatalError);
                }

            } // End if cell added cells empty
        } // End if face needs to be split into n
    } // End loop over all faces


    // PART 6: Adding faces

    // 6.1. Existing faces on special patches (empty or wedge) that get split
    //      (into n faces where n is the number of points or edges)
    // 6.2. Existing faces not on special patches that get split into two
    // 6.3. Existing faces that do not get split but only edges get split
    // 6.4. Existing faces that do not get split but get new owner/neighbour
    // 6.5. New internal faces inside split cells

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Marking faces to be handled"
            << endl;
    }

    // Get all faces to split:
    // a) All faces of a cell being split
    // b) All faces on special patch that are being split
    // c) All faces not on special patch that are being split
    // d) Both faces of an edge that is being split
    // Note: although a bit redundant, loop over everything above
    boolList facesToSplit(mesh_.nFaces(), false);

    // Also collect all split faces which will be needed in 6.3
    boolList allSplitFaces(mesh_.nFaces(), false);

    // Get edge faces
    const labelListList& meshEdgeFaces = mesh_.edgeFaces();

    // a) All faces of a cell that is being split
    forAll(refineCellsMask, cellI)
    {
        if (refineCellsMask.get(cellI))
        {
            const cell& cFaces = meshCells[cellI];

            forAll(cFaces, i)
            {
                facesToSplit[cFaces[i]] = true;
            }
        }
    }

    // b) All faces on special patch that are being split
    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] > -1)
        {
            // Mark face in both lists
            facesToSplit[faceI] = true; // Used through 6.1-6.5
            allSplitFaces[faceI] = true; // Used in 6.3
        }
    }

    // c) All faces not on special patch that are being split
    forAll(splitFacesIntoTwo, i)
    {
        // Get face index
        const label& faceI = splitFacesIntoTwo[i];

        // Mark face in both lists
        facesToSplit[faceI] = true;
        allSplitFaces[faceI] = true;
    }

    // d) Both faces of an edge that is being split
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


    // PART 6.1. Add/modify faces for each face on special patch that is being
    // split

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Splitting faces on special patches (empty or wedge)" << endl;
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
            forAll(f, fpI)
            {
                const label& pointI = f[fpI];

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
                        fpI,
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
                        fpI,
                        faceVerts
                    );

                    // Transfer dynamic list to a face (ordinary list)
                    newFace.transfer(faceVerts);
                    faceVerts.clear();

                    // Set new owner/neighbour indices based on split cells
                    label own, nei;
                    setNewFaceNeighbours
                    (
                        pointCellToAddedCellMap,
                        cellAddedCells,
                        faceI,
                        pointI, // Anchor point index

                        own,
                        nei
                    );

                    if (debug)
                    {
                        // Check orientation of the split face for debugging
                        checkNewFaceOrientation(meshMod, faceI, newFace);
                    }


                    // Finally insert the modification/addition instruction into
                    // the topo changer engine
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

        } // End if this face needs to be split
    } // End for all faces


    // PART 6.2. Add/modify faces for each face not on special patch that is
    // being split into two

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Splitting faces not on special patches (empty or wedge)"
            << endl;
    }

    // Loop through faces that are not on special patch. These will be split
    // into two faces only
    forAll (splitFacesIntoTwo, i)
    {
        // Get face index
        const label& faceI = splitFacesIntoTwo[i];

        // Check whether the face is marked for splitting. A bit redundant but
        // will be left for clarity
        if (facesToSplit[faceI])
        {
            // Face has not been split, grab the face
            const face& f = meshFaces[faceI];

            // Additional sanity check
            if (f.size() != 4)
            {
                FatalErrorInFunction
                    << "The original face has: " << f.size() << " points,"
                    << " while it should have exactly 4 points in order"
                    << " to split it in two."
                    << " faceI: " << faceI
                    << abort(FatalError);
            }

            // Flag to control whether the original faceI has been used
            // Note: original face gets modified, other face gets added
            bool modifiedFace = false;

            // Get anchor level for the face
            const label& anchorLevel = faceAnchorLevel[faceI];

            // Mark visited points to avoid adding the face twice
            boolList visitedPoint(f.size(), false);

            // New face always has four points/edges for 2D face splitting of
            // a face that is not on special patch
            face newFace(4);

            // Loop through all points of original face
            forAll (f, fpI)
            {
                // Get point index
                const label& pointI = f[fpI];

                if
                (
                    !visitedPoint[fpI]
                 && pointLevel_[pointI] <= anchorLevel
                )
                {
                    // This point is anchor and it hasn't been visited yet,
                    // start collecting face

                    // Collect the new face in the following order:
                    // 1. This point
                    // 2. Edge mid points for the edge that contains this point
                    //    and the edge that is on the other side
                    // 3. Other point (on the other side)

                    // 1. Set this point and mark it as visited
                    newFace[0] = pointI;
                    visitedPoint[fpI] = true;

                    // 2. Get edge mid point for edge containing this point
                    // Fetch the two edges on both sides
                    const Pair<label>& edgesOnOppositeSides =
                        splitFacesEmptyEdges[i];

                    // Get the edge indices and edges
                    const label& edgeIndexI = edgesOnOppositeSides.first();
                    const label& edgeIndexJ = edgesOnOppositeSides.second();

                    const edge& edgeI = meshEdges[edgeIndexI];
                    const edge& edgeJ = meshEdges[edgeIndexJ];

                    // Additional sanity check
                    if
                    (
                        (edgeMidPoint[edgeIndexI] == -1)
                     || (edgeMidPoint[edgeIndexJ] == -1)
                    )
                    {
                        // Edges are not marked for refinement, issue an error
                        FatalErrorInFunction
                            << "Trying to split a face into two, but"
                            << " edges on special patches (empty or wedge)"
                            << " are not properly set."
                            << nl
                            << "Edge: " << edgeIndexI << " with new point: "
                            << edgeMidPoint[edgeIndexI]
                            << " " << mesh_.points()[edgeI.start()]
                            << " " << mesh_.points()[edgeI.end()]
                            << nl
                            << "Edge: " << edgeIndexJ << " with new point: "
                            << edgeMidPoint[edgeIndexJ]
                            << " " << mesh_.points()[edgeJ.start()]
                            << " " << mesh_.points()[edgeJ.end()]
                            << abort(FatalError);
                    }

                    if ((edgeI.start() == pointI) || (edgeI.end() == pointI))
                    {
                        // Current point is on edgeI, set edgeI midpoint and
                        // then edgeJ midpoint
                        newFace[1] = edgeMidPoint[edgeIndexI];
                        newFace[2] = edgeMidPoint[edgeIndexJ];
                    }
                    else if
                    (
                        (edgeJ.start() == pointI) || (edgeJ.end() == pointI)
                    )
                    {
                        // Current is on edgeJ, set edgeJ midpoint and then
                        // edgeI midpoint
                        newFace[1] = edgeMidPoint[edgeIndexJ];
                        newFace[2] = edgeMidPoint[edgeIndexI];
                    }
                    else
                    {
                        // Point not on either of edges, issue an error
                        FatalErrorInFunction
                            << "Trying to split a face into two, but"
                            << " the point: " << pointI << " can't be found"
                            << " on either of edges. "
                            << nl
                            << "Edge: " << edgeIndexI << " with new point: "
                            << edgeMidPoint[edgeIndexI]
                            << nl
                            << "Edge: " << edgeIndexJ << " with new point: "
                            << edgeMidPoint[edgeIndexJ]
                            << abort(FatalError);
                    }

                    // At this point, we have added three points: original
                    // point, first edge mid point and second edge mid point.

                    // 3. Add the other point
                    // Get point edges for this point
                    const labelList& pEdges = meshPointEdges[pointI];

                    // Loop through all edges
                    forAll (pEdges, peI)
                    {
                        // Get the edge index
                        const label& pointEdgeI = pEdges[peI];

                        if (!edgeOnPatchToCut[pointEdgeI])
                        {
                            // Edge is not on special patch (empty or wedge),
                            // therefore this is the edge we're looking for.
                            // Collect the other point
                            const label pointJ =
                                meshEdges[pointEdgeI].otherVertex(pointI);

                            // Insert the point into the face at the last
                            // location
                            newFace[3] = pointJ;

                            // Mask local point index as visited by going
                            // through the face again
                            forAll (f, fpJ)
                            {
                                if (f[fpJ] == pointJ)
                                {
                                    // Found local index of the point, mask it
                                    visitedPoint[fpJ] = true;
                                }
                            }
                        }
                    }

                    // The face is now complete, set new owner/neighbour indices
                    // based on split cells
                    label own, nei;

                    // Set new face owner/neighbour pair
                    setNewFaceNeighbours
                    (
                        pointCellToAddedCellMap,
                        cellAddedCells,
                        faceI,
                        pointI, // Anchor point index

                        own,
                        nei
                    );

                    // We need to revert the face if the edge between this point
                    // and the next point is not split. This follows from
                    // definition of face as ordered set of points (defining
                    // orientation) and the splitting procedure. Note: edge
                    // ordering in face is the same as point ordering so the
                    // point index can be used as first face edge index
                    if (edgeMidPoint[meshFaceEdges[faceI][fpI]] == -1)
                    {
                        newFace = newFace.reverseFace();
                    }

                    if (debug)
                    {
                        // Check orientation of the split face for debugging
                        checkNewFaceOrientation(meshMod, faceI, newFace);
                    }


                    // Finally insert the modification/addition instruction into
                    // the topo changer engine
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
                } // End if point is anchored and has not been visited
            } // End loop over all face points

            // Mark face as handled
            facesToSplit[faceI] = false;

        } // End if face is split
    } // End for all faces that should be split into two


    // PART 6.3. Modify faces that do not get split but have edges that are
    // being split

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Modifying faces with split edges"
            << endl;
    }

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

                // Check whether this is not a face that's been split and that
                // the face has not been handled yet. The second check is
                // necessary since we go through edge faces instead of just
                // faces
                if (!allSplitFaces[faceI] && facesToSplit[faceI])
                {
                    // This is unsplit face that has not been handled

                    // Get face and face edges
                    const face& f = meshFaces[faceI];
                    const labelList& fEdges = meshFaceEdges[faceI];

                    // Create a dynamic list containing new face vertices
                    DynamicList<label> newFaceVerts(f.size());

                    // Append all original points and all edge mid points
                    forAll(f, fpI)
                    {
                        newFaceVerts.append(f[fpI]);

                        const label edgeI = fEdges[fpI];

                        if (edgeMidPoint[edgeI] > -1)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    // Create a face from dynamic list by transfer
                    face newFace(move(newFaceVerts));


                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    const label anchorFpI = findMinLevel(f);

                    label own, nei;
                    setNewFaceNeighbours
                    (
                        pointCellToAddedCellMap,
                        cellAddedCells,
                        faceI,
                        f[anchorFpI], // Anchor point index

                        own,
                        nei
                    );


                    if (debug)
                    {
                        // Check orientation of the new face for debugging
                        checkNewFaceOrientation(meshMod, faceI, newFace);
                    }

                    // Modify the face
                    modifyFace(meshMod, faceI, newFace, own, nei);

                    // Mark face as handled
                    facesToSplit[faceI] = false;

                } // End if unsplit, unhandled face
            } // End for all edge faces
        } // End if edge has been cut
    } // End for all edges


    // PART 6.4: Modify faces that do not get split but whose owner/neighbour
    // change due to splitting

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Changing owner/neighbour for otherwise unaffected faces"
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
            label anchorFpI = findMinLevel(f);

            label own, nei;
            setNewFaceNeighbours
            (
                pointCellToAddedCellMap,
                cellAddedCells,
                faceI,
                f[anchorFpI], // Anchor point

                own,
                nei
            );

            // Modify the face, changing owner and neighbour
            modifyFace(meshMod, faceI, f, own, nei);

            // Mark face as handled
            facesToSplit[faceI] = false;

        } // End if the face needs to be handled
    } // End for all faces


    // PART 6.5. Add new internal faces inside split cells

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Adding new internal faces for split cells"
            << endl;
    }

    // Mark-up filed for visited cells (since we are going through faces)
    boolList cellsToSplit(mesh_.nCells(), true);

    // Loop through faces in the same way as we did when we were adding
    // cells. This order is important since it ensures easy determination of
    // owner/neighbour cells for new faces
    forAll(faceMidPoint, faceI)
    {
        // Get owner of the face. For face on special patch (empty or wedge)
        const label& cellI = meshFaceOwner[faceI];

        // Check whether this face has been split and whether the cell has been
        // handled (internal faces already created for this cell)
        if (faceMidPoint[faceI] > -1 && cellsToSplit[cellI])
        {
            // Face is split and hasn't been visited yet. Get the face and edges
            const face& f = meshFaces[faceI];
            const labelList& fEdges = meshFaceEdges[faceI];

            // Get anchor points for this cell and cell added cells
            const labelList& cAnchors = cellAnchorPoints[cellI];
            const labelList& cAdded = cellAddedCells[cellI];

            // Count number of added faces (helper variable to determine
            // owner/neighbour)
            label nAddedFaces = 0;

            // Loop through face points
            forAll (f, fpI)
            {
                // Get point index
                const label& pointI = f[fpI];

                // If this point is not an anchor, it has already been handled
                // (by going through anchors), skip
                if (findIndex(cAnchors, pointI) == -1)
                {
                    continue;
                }

                // Get corresponding edge index and edge (between fpI and
                // fpI + 1) by definition of faceEdges
                const label& edgeI = fEdges[fpI];
                const edge& e = meshEdges[edgeI];

                // Grab other point
                const label& pointJ = e.otherVertex(pointI);

                if (pointJ == -1)
                {
                    // If pointJ is equal to -1, this means that the pointI
                    // was not found on edge, something went wrong
                    FatalErrorInFunction
                        << "Point: " << pointI << " not found on edge: "
                        << edgeI << nl
                        << "Looping through face points and face edges did"
                        << " not ensure synchronous behaviour."
                        << abort(FatalError);
                }

                // Create the new face
                face newFace(4);

                // Note: there are three possible variants:
                //   i) Edge is split and the other point is anchor. Collection
                //      of the face starts from edgeMidPoint
                //  ii) Edge is split and the other point is not an
                //      anchor. Collection of the face starts from other point
                // iii) Edge is not split and the other point is not an
                //      anchor. Collection of the face starts from other point
                // Variants ii) and iii) can be handled together, while variant
                // i) has to be handled separately.

                // Whether the edge is split
                const bool isEdgeSplit = edgeMidPoint[edgeI] > -1;

                // Whether the other point is anchor or not
                const bool isOtherEdgePointAnchor
                    = findIndex(cAnchors, pointJ) > -1;

                // Check if the edge is split and whether the other edge point
                // is an anchor
                if (isEdgeSplit && isOtherEdgePointAnchor)
                {
                    // Variant i) Edge is split and other edge point is anchor

                    // Create the new face and start collecting points
                    // a) edgeMidPoint for this edge
                    // b) faceMidPoint for this face
                    // c) faceMidPoint for the face on the other side
                    // d) edgeMidPoint for the edge on the other side

                    // a) edgeMidPoint for this edge
                    newFace[0] = edgeMidPoint[edgeI];

                    // b) and c): adding both face mids
                    addFaceMids
                    (
                        faceMidPoint,
                        faceOnPatchToCut,
                        faceI,
                        cellI,
                        newFace
                    );

                    // d) edgeMidPoint for the edge on the other side
                    // The other edge is uniquely defined as the edge on special
                    // patch (empty or wedge) sharing the same face as this edge

                    // Get the edge faces
                    const labelList& eFaces = meshEdgeFaces[edgeI];

                    // Loop through edge faces
                    forAll(eFaces, i)
                    {
                        // Get the face and check whether it is on special patch
                        const label& faceK = eFaces[i];

                        if (!faceOnPatchToCut[faceK])
                        {
                            // Found the face, need to search its edges
                            const labelList& otherFaceEdges =
                                meshFaceEdges[faceK];

                            forAll(otherFaceEdges, j)
                            {
                                // Get the edge
                                const label& edgeJ = otherFaceEdges[j];

                                if (edgeOnPatchToCut[edgeJ] && (edgeI != edgeJ))
                                {
                                    // Edge is on special patch (empty or
                                    // wedge), this must be the one we are
                                    // looking for. Add its midpoint and double
                                    // check if it is valid
                                    if (edgeMidPoint[edgeJ] > -1)
                                    {
                                        newFace[3] = edgeMidPoint[edgeJ];
                                        break;
                                    }
                                    else
                                    {
                                        FatalErrorInFunction
                                            << "Other edge: "
                                            << edgeJ
                                            << " has not been selected for"
                                            << " splitting, while the edge on"
                                            << " original side: "
                                            << edgeI
                                            << " has been selected."
                                            << abort(FatalError);
                                    }
                                } // End if this is our "other" edge
                            } // End for all other (non special patch: empty or
                              // wedge) face edges

                            // Break out since we must have found the candidate
                            break;

                        } // End if face not on special patch (empty or wedge)

                    } // End for all edge faces

                } // End if this edge is split and the other point is anchor
                else if (!isOtherEdgePointAnchor)
                {
                    // Variants ii) and iii). Either the edge is split and the
                    // other point is not an anchor or the edge is not split and
                    // the other point is not an anchor

                    // Create the new face and start collecting points
                    // a) other point of this edge
                    // b) faceMidPoint for this face
                    // c) faceMidPoint for the face on the other side
                    // d) other point on the other side

                    // a) other point of this edge
                    newFace[0] = pointJ;

                    // b) and c): adding both face mids
                    addFaceMids
                    (
                        faceMidPoint,
                        faceOnPatchToCut,
                        faceI,
                        cellI,
                        newFace
                    );

                    // d) other point on the other side
                    // The other point is uniquely defined as the other point of
                    // the edge of this point which is not on special patch

                    // Get point edges
                    const labelList& pEdges = meshPointEdges[pointJ];

                    // Loop through all edges
                    forAll(pEdges, i)
                    {
                        // Get the edge index
                        const label& edgeJ = pEdges[i];

                        if (!edgeOnPatchToCut[edgeJ])
                        {
                            // Found our edge, set the point on the other side
                            // of the edge as the last point in face
                            newFace[3] = meshEdges[edgeJ].otherVertex(pointJ);
                            break;
                        }
                    } // End loop over all point edges

                } // End if the other point is not an anchor
                else
                {
                    // The edge is not split and the other point is an
                    // anchor. This should never happen
                    FatalErrorInFunction
                        << "Attempted to create internal face for an edge that"
                        << " is not split and the other point that is an anchor."
                        << nl
                        << "Cell: " << cellI
                        << ", point: " << pointI
                        << ", other edge point: " << pointJ
                        << nl
                        << "Anchor points for cell are: " << cAnchors
                        << abort(FatalError);
                }

                // Now we have the face defined, set owner and neighbour.
                // Note: owner and neighbour are uniquely defined since we have
                // gone through the face in the same way as we did while adding
                // cells. This ensured easy definition of owner/neighbour cells
                label own = cAdded[nAddedFaces];
                label nei = cAdded[cAdded.fcIndex(nAddedFaces)];

                // According to the definition of adding faces, the first n - 1
                // faces need to be reverted, while the last one is correctly
                // oriented

                // ***** THIS LINE SEEMS TO CAUSE PROBLEMS ******
                if (nAddedFaces < cAdded.size() - 1)
                {
                    newFace = newFace.reverseFace();
                }
                else
                {
                    Swap(own, nei);
                }


                // Debug: check orientation
                if (debug)
                {
                    // Get owner/neighbour points
                    point ownPt, neiPt;

                    if (nAddedFaces < cAdded.size() - 1)
                    {
                        // Original owner/neighbour
                        ownPt = meshPoints[pointI];
                        neiPt = meshPoints[pointJ];
                    }
                    else
                    {
                        // Flipped owner/neighbour for last face
                        ownPt = meshPoints[pointJ];
                        neiPt = meshPoints[pointI];
                    }
                    checkInternalOrientation
                    (
                        meshMod,
                        cellI,
                        faceI,
                        ownPt,
                        neiPt,
                        newFace
                    );
                }

                // Finally, add the face. Note: ignoring return of new face
                // index from meshMod.setAction(polyAddFace(...)) call
                addInternalFace
                (
                    meshMod,
                    faceI,
                    nAddedFaces < cAdded.size() - 1 ? pointI : pointJ,
                    newFace,
                    own,
                    nei
                );

                // Increment number of added faces
                nAddedFaces++;

            } // End loop over all point (and edges) of the face

            // Finished adding internal faces. Mark the cell as handled
            cellsToSplit[cellI] = false;

        } // End if face is split into n and cell has not been handled
    } // End for all faces

    // Debug: check minimum point index of added points, needs to be equal to
    // number of points in the original mesh
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

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


void Foam::prismatic2DRefinement::setUnrefinement
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

    {
        labelList newCellLevel(cellLevel_);
        // Update refinementLevelIndicator for all cells that will be unrefined
        forAll(splitPointsToUnrefine, i)
        {
            // Get point cells and mark them for unrefinement
            const labelList& pCells = meshPointCells[splitPointsToUnrefine[i]];

            forAll(pCells, j)
            {
                newCellLevel[pCells[j]] = cellLevel_[pCells[j]] - 1;
            }
        }
        cellLevel_.transfer(newCellLevel);
    }

    // Create lists needed by face remover
    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    // Memory management
    {
        // Mark faces on special patches (empty or wedge) to exclude them
        boolList faceOnPatchToCut(mesh_.nFaces(), false);

        // Get boundary
        const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();

        // Loop through all patches
        forAll (boundaryMesh, patchI)
        {
            // Get current patch
            const polyPatch& curPatch = boundaryMesh[patchI];

            // Check whether this patch is special (empty or wedge)
            if (isA<emptyPolyPatch>(curPatch) || isA<wedgePolyPatch>(curPatch))
            {
                // Get start and end face labels
                const label startFaceI = curPatch.start();
                const label endFaceI = startFaceI + curPatch.size();

                // Mark all the faces and edges on the patch
                for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
                {
                    // Mark face
                    faceOnPatchToCut[faceI] = true;
                }
            }
        }


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
                // Get face index
                const label& faceI = pFaces[j];

                if (!faceOnPatchToCut[faceI])
                {
                    // Face is not on special patch, insert it into hash set
                    splitFaces.insert(faceI);
                }
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prismatic2DRefinement::prismatic2DRefinement
(
    const polyMesh& mesh,
    const dictionary& dict,
    const bool read
)
:
    refinement(mesh, dict, read)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::prismatic2DRefinement::~prismatic2DRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::labelList Foam::prismatic2DRefinement::consistentUnrefinement
(
    const labelList& unrefinementPointCandidates,
    const bool maxSet
) const
{
    if (debug)
    {
        InfoInFunction
            << "Setting split points to unrefine." << endl;
    }

    // Get necessary mesh data
    const label nPoints = mesh_.nPoints();

    // PART 1: Mark all split points in the mesh (points that can be unrefined)
    PackedBoolList splitPointsMarkup(nPoints);

    // Algorithm: split point is uniquely defined as a point that:
    // 1. Has pointLevel_ > 0 (obviously),
    // 2. A point that has the same pointLevel_ as ALL of the points of its
    //    edges. In other words, for each point, we will look through all the
    //    edges of the point. For each edge, we will visit both points and
    //    check point levels. All point levels must be the same for this point
    //    candidate to be a split point. This is quite useful since there is no
    //    need to store the refinement history

    // Get necessary mesh data
    const edgeList& meshEdges = mesh_.edges();
    const labelListList& meshPointEdges = mesh_.pointEdges();

    // Loop through all points
    forAll (meshPointEdges, pointI)
    {
        // Get point level of this point
        const label& centralPointLevel = pointLevel_[pointI];

        if (!centralPointLevel)
        {
            // Point can't be unrefined as its level is either 0 or
            // invalid. Continue immediately
            continue;
        }

        // Flag to see whether this is a split point candidate
        bool splitPointCandidate = true;

        // Get all edge labels for this point
        const labelList& pEdges = meshPointEdges[pointI];

        // Loop through all point edges
        forAll (pEdges, i)
        {
            // Get edge index and the edge
            const label& edgeI = pEdges[i];
            const edge& curEdge = meshEdges[edgeI];

            // Loop through both points of the edge
            forAll (curEdge, j)
            {
                // Get point index
                const label& pointJ = curEdge[j];

                if (pointLevel_[pointJ] != centralPointLevel)
                {
                    // Point levels are different, this can't be a split point,
                    // set flag to false and break immediatelly
                    splitPointCandidate = false;
                    break;
                }
                // else: this is still potential split point candidate so
                //       there's nothing to do
            } // End for both points of this edge

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
            splitPointsMarkup.set(pointI);
        }
    }

    // Note: if there is no dynamic load balancing, points at the processor
    // boundary cannot be split points by definition. However, in dynamic load
    // balancing runs, it is possible that a split point end on processor
    // boundary, in which case we will simply avoid (actually delay) unrefining
    // until this becomes internal point again. VV, 4/Jul/2018.

    // Get boundary mesh and mesh faces
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const faceList& meshFaces = mesh_.faces();

    // Loop through all patches
    forAll (bMesh, patchI)
    {
        const polyPatch& patch = bMesh[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            // Get patch start
            const label startIndex = patch.start();

            // Loop through all the faces
            forAll (patch, i)
            {
                // Get global face index and face
                const label faceI = startIndex + i;
                const face& f = meshFaces[faceI];

                // Make sure that we don't split around point at all points of
                // the processor patch faces
                forAll (f, fpI)
                {
                    splitPointsMarkup.unset(f[fpI]);
                }
            }
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

        if (splitPointsMarkup.get(pointI))
        {
            // This is a split point, mark it for unrefinement
            splitPointsToUnrefine.set(pointI);
        }
    }


    // PART 3: Ensure face consistent (2:1 constraint) and possibly edge
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
        PackedBoolList cellsToUnrefine(nCells);

        // Loop through all split points to unrefine
        forAll (splitPointsToUnrefine, pointI)
        {
            if (splitPointsToUnrefine[pointI])
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

        // Check for 2:1 face based consistent unrefinement. Updates
        // cellsToUnrefine and returns number of removed cells from unrefinement
        // in this iteration
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
                    if (!cellsToUnrefine.get(pCells[i]))
                    {
                        // Cell must not be refined, remove point from
                        // unrefinement as well
                        splitPointsToUnrefine.unset(pointI);
                        break;
                    }
                }
            }
        }
    }

    // Convert back to labelList.
    label nSet = splitPointsToUnrefine.count();
    labelList newPointsToUnrefine(nSet);
    nSet = 0;

    forAll(splitPointsToUnrefine, pointi)
    {
        if (splitPointsToUnrefine.get(pointi))
        {
            newPointsToUnrefine[nSet++] = pointi;
        }
    }

    return newPointsToUnrefine;
}


// ************************************************************************* //
