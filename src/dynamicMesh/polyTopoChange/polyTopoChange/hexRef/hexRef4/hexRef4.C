/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "hexRef4.H"

#include "emptyPolyPatch.H"

#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef4, 0);
    addToRunTimeSelectionTable(hexRef, hexRef4, mesh);
    addToRunTimeSelectionTable(hexRef, hexRef4, levelsHist);
    addToRunTimeSelectionTable(hexRef, hexRef4, levels);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Check whether pointi is an anchor on celli.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRef4::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label celli,
    const label facei,
    const label pointi
) const
{
    if (cellAnchorPoints[celli].size())
    {
        label index = findIndex(cellAnchorPoints[celli], pointi);

        if (index != -1)
        {
            if (index >= 4) //AB....
            {
                if (index == 4)
                {
                    index = 8;
                }
                index = 8 - index;
            } //AB
            return cellAddedCells[celli][index];
        }


        // pointi is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[facei];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[celli], f[fp]);

            if (index != -1)
            {
                if (index >= 4) //AB....
                {
                if (index == 4)
                    {
                        index = 8;
                    }
                    index = 8 - index;
                } //AB
                return cellAddedCells[celli][index];
            }
        }

        // Problem.
        dumpCell(celli);
        Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
            << endl;

        FatalErrorInFunction
            << "Could not find point " << pointi
            << " in the anchorPoints for cell " << celli << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return celli;
    }
}


// Internal faces are one per edge between anchor points. So one per midPoint
// between the anchor points. Here we store the information on the midPoint
// and if we have enough information:
// - two anchors
// - two face mid points
// we add the face. Note that this routine can get called anywhere from
// two times (two unrefined faces) to four times (two refined faces) so
// the first call that adds the information creates the face.
Foam::label Foam::hexRef4::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& edgeMidPoint,
    const label celli,
    const label facei,
    const bool faceOrder,
    const label edgeMidPointi,
    const label anchorPointi,
    const label faceMidPointi,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    // See if need to store anchors.

    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointi);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointi, edge(anchorPointi, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointi != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointi;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }


    // Check if this call of storeMidPointinfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors)
    {
        const cell& cFaces = mesh_.cells()[celli];
        label face1 = -1;

        forAll(cFaces, i)
        {
            label faceJ = cFaces[i];

            if (faceMidPoint[faceJ] != faceMidPointi
                && faceMidPoint[faceJ] >= 0
                && faceMidPoint[faceJ] != 1234567890) // could be replace by faceJ == type "empty"
            {
                face1 = faceJ;  // needed to take the other face points
            }
        }

        const edge& anchors = midPointToAnchors[edgeMidPointi];
        label index = findIndex(cellAnchorPoints[celli], anchorPointi);

        if (index == 0)
        {
            index = 4;
        }
        if (index == 4)
        {
            index = 8;
        }

        label point1 = cellAnchorPoints[celli][8 - index];
        label edgeMidPointj = -1;

        const face& f = mesh_.faces()[face1];       //other face
        const labelList& fEdges = mesh_.faceEdges(face1);

        DynamicList<label> newFaceVerts(4);

        if (faceOrder == (mesh_.faceOwner()[facei] == celli))
        {
            label anch = findIndex(f, point1);

            if (pointLevel_[f[f.rcIndex(anch)]] <= cellLevel_[celli])
            {
                label edgeJ = fEdges[f.rcIndex(anch)];
                edgeMidPointj = edgeMidPoint[edgeJ];
            }
            else
            {
                label edgeMid = findLevel(face1, f, f.rcIndex(anch), false, cellLevel_[celli] +1);
                edgeMidPointj = f[edgeMid];
            }

            newFaceVerts.append(faceMidPointi);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointi,  // edge between faceMid and
                edgeMidPointi,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                edgeMidPointj,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointj);
            newFaceVerts.append(faceMidPoint[face1]);
        }
        else
        {
            label anch = findIndex(f, point1);

            if (pointLevel_[f[f.fcIndex(anch)]] <= cellLevel_[celli])
            {
                label edgeJ = fEdges[anch];
                edgeMidPointj = edgeMidPoint[edgeJ];
            }
            else
            {
                label edgeMid = findLevel
                (
                    face1,
                    f,
                    f.fcIndex(anch),
                    false,
                    cellLevel_[celli] +1
                );
                edgeMidPointj = f[edgeMid];
            }

            newFaceVerts.append(edgeMidPointj);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointj,
                edgeMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                faceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointi);
            newFaceVerts.append(faceMidPoint[face1]);
        }

        face newFace;
        newFace.transfer(newFaceVerts);

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            celli,
            facei,
            anchorPointi
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            celli,
            facei,
            anchors.otherVertex(anchorPointi)
        );


        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = mesh_.points()[anchorPointi];
            neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace.flip();

            ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
            neiPt = mesh_.points()[anchorPointi];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = mesh_.points()[anchorPointi];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
            }
            else
            {
                ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
                neiPt = mesh_.points()[anchorPointi];
            }

            checkInternalOrientation
            (
                meshMod,
                celli,
                facei,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            meshMod,
            facei,
            anchorPointi,
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


// Creates all the 4 internal faces for celli.
void Foam::hexRef4::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label celli,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[celli];
    const label cLevel = cellLevel_[celli];

    // From edge mid to anchor points
    Map<edge> midPointToAnchors(8);
    // From edge mid to face mids
    Map<edge> midPointToFaceMids(8);

    // Storage for on-the-fly addressing
    DynamicList<label> storage;


    // Running count of number of internal faces added so far.
    label nFacesAdded = 0;

    forAll(cFaces, i)
    {
        label facei = cFaces[i];

        const face& f = mesh_.faces()[facei];
        const labelList& fEdges = mesh_.faceEdges(facei, storage);

        // We are on the celli side of face f. The face will have 1 or 4
        // cLevel points and lots of higher numbered ones.

        label faceMidPointi = -1;

        label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.
            // In 2D case this should never happen, as the shared faces with uncutted mesh
            // can be only the ones splitted in two which have two anchors and two middle edge points.

            Info << "Should never happen: nAnchors == 1" << endl;

            // Find the one anchor.
            label anchorFp = -1;

            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel+1 point:
            // these two lines just repeat the same steps two times to reach the second cLevel+1
            label edgeMid = findLevel
            (
                facei,
                f,
                f.fcIndex(anchorFp),
                true,
                cLevel+1
            );

            // to be checked!!!
            label faceMid = findLevel
            (
                facei,
                f,
                f.fcIndex(edgeMid),
                true,
                cLevel+1
            );

            faceMidPointi = f[faceMid];
        }
        else if (nAnchors == 2)
        {
            // Only two anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.

            faceMidPointi = 1234567890;        // again not nice

        }
        else if (nAnchors == 4)
        {
            // There is no face middle yet but the face will be marked for
            // splitting. For non-empty faces, faceMidPoint[facei] will be 1234567890,
            // which is not a point ID. This is checked in storeMidPointinfo.

            faceMidPointi = faceMidPoint[facei];
        }
        else
        {
            dumpCell(mesh_.faceOwner()[facei]);
            if (mesh_.isInternalFace(facei))
            {
                dumpCell(mesh_.faceNeighbour()[facei]);
            }

            FatalErrorInFunction
                << "nAnchors:" << nAnchors
                << " facei:" << facei
                << abort(FatalError);
        }



        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointinfo will collect
        // all the info and combine it all.
        // Only empty faces are considered as they are split like the 3D case.

        if(faceMidPoint[facei] != 1234567890 && faceMidPointi != 1234567890)
        {
            forAll(f, fp0)
            {
                label point0 = f[fp0];

                if (pointLevel_[point0] <= cLevel)
                {
                    // Anchor.

                    // Walk forward
                    // ~~~~~~~~~~~~
                    // to cLevel+1 or edgeMidPoint of this level.


                    label edgeMidPointi = -1;

                    label fp1 = f.fcIndex(fp0);

                    if (pointLevel_[f[fp1]] <= cLevel)
                    {
                        // Anchor. Edge will be split.
                        label edgeI = fEdges[fp0];

                        edgeMidPointi = edgeMidPoint[edgeI];

                        if (edgeMidPointi == -1)
                        {
                            dumpCell(celli);

                            const labelList& cPoints = mesh_.cellPoints(celli);

                            FatalErrorInFunction
                                << "cell:" << celli << " cLevel:" << cLevel
                                << " cell points:" << cPoints
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, cPoints)()
                                << " face:" << facei
                                << " f:" << f
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, f)()
                                << " faceAnchorLevel:" << faceAnchorLevel[facei]
                                << " faceMidPoint:" << faceMidPoint[facei]
                                << " faceMidPointi:" << faceMidPointi
                                << " fp:" << fp0
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // Search forward in face to clevel+1
                        // In 2D this should never be used as this case is verified with
                        // shared faces between refined and not-refined cells which have
                        // to be refined in this turn.

                        label edgeMid = findLevel(facei, f, fp1, true, cLevel+1);

                        edgeMidPointi = f[edgeMid];
                    }

                    label newFacei = storeMidPointInfo
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        cellMidPoint,
                        faceMidPoint,
                        edgeMidPoint,

                        celli,
                        facei,
                        true,                   // mid point after anchor
                        edgeMidPointi,          // edgemid
                        point0,                 // anchor
                        faceMidPointi,

                        midPointToAnchors,
                        midPointToFaceMids,
                        meshMod
                    );

                    if (newFacei != -1)
                    {
                        nFacesAdded++;

                        if (nFacesAdded == 4)
                        {
                            break;
                        }
                    }



                    // Walk backward
                    // ~~~~~~~~~~~~~

                    label fpMin1 = f.rcIndex(fp0);

                    if (pointLevel_[f[fpMin1]] <= cLevel)
                    {
                        // Anchor. Edge will be split.
                        label edgeI = fEdges[fpMin1];

                        edgeMidPointi = edgeMidPoint[edgeI];

                        if (edgeMidPointi == -1)
                        {
                            dumpCell(celli);

                            const labelList& cPoints = mesh_.cellPoints(celli);

                            FatalErrorInFunction
                                << "cell:" << celli << " cLevel:" << cLevel
                                << " cell points:" << cPoints
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, cPoints)()
                                << " face:" << facei
                                << " f:" << f
                                << " pointLevel:"
                                << UIndirectList<label>(pointLevel_, f)()
                                << " faceAnchorLevel:" << faceAnchorLevel[facei]
                                << " faceMidPoint:" << faceMidPoint[facei]
                                << " faceMidPointi:" << faceMidPointi
                                << " fp:" << fp0
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // Search back to clevel+1
                        // In 2D this should never be used as this case is verified with
                        // shared faces between refined and not-refined cells which have
                        // to be refined in this turn.

                        label edgeMid = findLevel
                        (
                            facei,
                            f,
                            fpMin1,
                            false,
                            cLevel+1
                        );

                        edgeMidPointi = f[edgeMid];
                    }

                    newFacei = storeMidPointInfo
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        cellMidPoint,
                        faceMidPoint,
                        edgeMidPoint,

                        celli,
                        facei,
                        false,                  // mid point before anchor
                        edgeMidPointi,          // edgemid
                        point0,                 // anchor
                        faceMidPointi,

                        midPointToAnchors,
                        midPointToFaceMids,
                        meshMod
                    );

                    if (newFacei != -1)
                    {
                        nFacesAdded++;

                        if (nFacesAdded == 4)
                        {
                            break;
                        }
                    }
                }   // done anchor
            } // done face
        }   // done check empty condition

        if (nFacesAdded == 4)
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef4::hexRef4(const polyMesh& mesh, const bool readHistory)
:
    hexRef(mesh, readHistory)
{}


// Construct from components
Foam::hexRef4::hexRef4
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory& history,
    const scalar level0Edge
)
:
    hexRef(mesh, cellLevel, pointLevel, history, level0Edge)
{}


// Construct from components
Foam::hexRef4::hexRef4
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    hexRef(mesh, cellLevel, pointLevel, level0Edge)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Top level driver to insert topo changes to do all refinement.
Foam::labelListList Foam::hexRef4::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        //checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();


    // New point/cell level. Copy of pointLevel for existing points.
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
        Pout<< "hexRef4::setRefinement :"
            << " Allocating " << cellLabels.size() << " cell midpoints."
            << endl;
    }


    // Mid point per refined cell.
    // -1 : not refined
    // >=0: label of mid point.
    labelList cellMidPoint(mesh_.nCells(), -1);
    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];
        cellMidPoint[celli] = 1234567890;    // mark need for splitting
    }


    boolList isDivisibleFace(mesh_.nFaces(), false);
    boolList isDivisibleEdge(mesh_.nEdges(), false);

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const label & patchID = mesh_.boundaryMesh().whichPatch(facei);

//        if (isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchID])
        if (mesh_.boundaryMesh()[patchID].type() == "empty")
        {
            isDivisibleFace[facei] = true;
            const labelList& fEdges = mesh_.faceEdges(facei);

            forAll(fEdges, i)
            {
                label edgeJ = fEdges[i];
                isDivisibleEdge[edgeJ] = true;
            }
        }
    }


    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellLabels.size());

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                splitCells.insert(celli);
            }
        }

        Pout<< "hexRef4::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }



    // Split edges
    // ~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over cells to be refined or edges?
    forAll(cellMidPoint, celli)
    {
        if (cellMidPoint[celli] >= 0)
        {
            const labelList& cEdges = mesh_.cellEdges(celli);
            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];
                const edge& e = mesh_.edges()[edgeI];
                if
                (
                 isDivisibleEdge[edgeI]
                 && pointLevel_[e[0]] <= cellLevel_[celli]
                 && pointLevel_[e[1]] <= cellLevel_[celli]
                )
                {
                    edgeMidPoint[edgeI] = 12345;    // mark need for splitting
                }
            }
        }
    }

    // Synchronize edgeMidPoint across coupled patches. Take max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
        );


        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    )
                  + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexRef4::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Calculate face level
    // ~~~~~~~~~~~~~~~~~~~~
    // (after splitting)

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label facei = 0; facei < mesh_.nFaces(); facei++)
    {
        faceAnchorLevel[facei] = faceLevel(facei);
    }

    // -1  : no need to split face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);


    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceAnchorLevel[facei] >= 0)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei];
            label newNeiLevel = neiLevel + (cellMidPoint[nei] >= 0 ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[facei]
             || newNeiLevel > faceAnchorLevel[facei]
            )
            {
                faceMidPoint[facei] = 1234567890;    // mark to be split (not a real nice way)
            }
        }
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceAnchorLevel[facei] >= 0)
            {
                label own = mesh_.faceOwner()[facei];
                label ownLevel = cellLevel_[own];
                label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

                if
                (
                    newOwnLevel > faceAnchorLevel[facei]
                 || newNeiLevel[i] > faceAnchorLevel[facei]
                )
                {
                    faceMidPoint[facei] = 1234567890;    // mark to be split (not really nice way)
                }
            }
        }
    }


    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );



    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: determine mid points and sync. See comment for edgeMids
        // above
        pointField bFaceMids
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            point(-GREAT, -GREAT, -GREAT)
        );

        forAll(bFaceMids, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceMidPoint[facei] >= 0)
            {
                bFaceMids[i] = mesh_.faceCentres()[facei];
            }
        }
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
        );

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0 && isDivisibleFace[facei])
            {
                // Face marked to be split. Replace faceMidPoint with actual
                // point label.

                const face& f = mesh_.faces()[facei];
                faceMidPoint[facei] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            facei < mesh_.nInternalFaces()
                          ? mesh_.faceCentres()[facei]
                          : bFaceMids[facei-mesh_.nInternalFaces()]
                        ),                          // point
                        f[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                // Determine the level of the corner points and midpoint will
                // be one higher.
                newPointLevel(faceMidPoint[facei]) = faceAnchorLevel[facei]+1;
            }
        }
    }

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", cellLabels.size());

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                splitFaces.insert(facei);
            }
        }

        Pout<< "hexRef4::setRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - cellMidPoint >= 0 : cell needs to be split.
    // - faceMidPoint >= 0 : face needs to be split.
    // - edgeMidPoint >= 0 : edge needs to be split.
    // - isDivisibleFace true : face needs to be split in 4 else in 2.
    //   for the 2 split case, faceMidPoint is >= 0 but not initalized with a point.


    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());
    {
        labelList nAnchorPoints(mesh_.nCells(), 0);
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                cellAnchorPoints[celli].setSize(8);
            }
        }

       forAll(cellMidPoint, celli)
       {
            const cell& cFaces = mesh_.cells()[celli];
            forAll(cFaces, i)
            {
                label facei = cFaces[i];
                const face& f = mesh_.faces()[facei];
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if
                    (
                        cellMidPoint[celli] >= 0
                        && isDivisibleFace[facei]
                        && pointLevel_[pointi] <= cellLevel_[celli]
                    )
                    {
                        if (nAnchorPoints[celli] == 8)
                        {
                            dumpCell(celli);
                            FatalErrorInFunction
                                << "cell " << celli
                                << " of level " << cellLevel_[celli]
                                << " uses more than 8 points of equal or"
                                << " lower level" << nl
                                << "Points so far:" << cellAnchorPoints[celli]
                                << abort(FatalError);
                        }
                        cellAnchorPoints[celli][nAnchorPoints[celli]++]
                            = pointi;
                    }
                }
            }
        }


        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                if (nAnchorPoints[celli] != 8)
                {
                    const labelList cPoints(mesh_.cellPoints(celli));

                    FatalErrorInFunction
                        << "cell " << celli
                        << " of level " << cellLevel_[celli]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << IndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Adding cells (1 per anchorPoint)"
            << endl;
    }

    // Per cell the 3 added cells (+ original cell)
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellAnchorPoints, celli)
    {
        const labelList& cAnchors = cellAnchorPoints[celli];

        if (cAnchors.size() == 8)
        {
            labelList& cAdded = cellAddedCells[celli];
            cAdded.setSize(4);

            // Original cell at 0
            cAdded[0] = celli;

            // Update cell level
            newCellLevel[celli] = cellLevel_[celli]+1;


            for (label i = 1; i < 4; i++)
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                                 // master point
                        -1,                                 // master edge
                        -1,                                 // master face
                        celli,                              // master cell
                        mesh_.cellZones().whichZone(celli)  // zone for cell
                    )
                );

                newCellLevel(cAdded[i]) = cellLevel_[celli]+1;
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but only edges get split
    // 3. existing faces that do not get split but get new owner/neighbour
    // 4. new internal faces inside split cells.

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, i)
                {
                    affectedFace.set(cFaces[i]);
                }
            }
        }

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                affectedFace.set(facei);
            }
        }

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const labelList& eFaces = mesh_.edgeFaces(edgeI);

                forAll(eFaces, i)
                {
                    affectedFace.set(eFaces[i]);
                }
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, facei)
    {
        if (faceMidPoint[facei] >= 0 && affectedFace.get(facei))
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[facei];

            // Has original facei been used (three faces added, original gets
            // modified)
            bool modifiedFace = false;
            label anchorLevel = faceAnchorLevel[facei];

            if(isDivisibleFace[facei])
            {
                face newFace(4);
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (pointLevel_[pointi] <= anchorLevel)
                    {
                        // point is anchor. Start collecting face.
                        DynamicList<label> faceVerts(4);
                        faceVerts.append(pointi);

                        // Walk forward to mid point.
                        // - if next is +2 midpoint is +1
                        // - if next is +1 it is midpoint
                        // - if next is +0 there has to be edgeMidPoint

                        walkFaceToMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        faceVerts.append(faceMidPoint[facei]);
                        walkFaceFromMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        // Convert dynamiclist to face.
                        newFace.transfer(faceVerts);
                        // Get new owner/neighbour
                        label own, nei;
                        getFaceNeighbours
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            facei,
                            pointi,          // Anchor point
                            own,
                            nei
                        );


                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];
                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }
                        if (!modifiedFace)
                        {
                            modifiedFace = true;
                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }
            }
            else
            {
                face newFace(2);
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    label nextpointi = f[f.fcIndex(fp)];
                    label edgeI = meshTools::findEdge (mesh_, pointi, nextpointi);
                    if (edgeMidPoint[edgeI] >=0)
                    {
                        DynamicList<label> faceVerts(4);
                        label pointJ = f[f.rcIndex(fp)];
                        faceVerts.append(pointi);
                        walkFaceToMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        walkFaceFromMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            f.rcIndex(fp),
                            faceVerts
                        );

                        faceVerts.append(pointJ);
                        newFace.transfer(faceVerts);
                        label own, nei;
                        getFaceNeighbours
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            facei,
                            pointi,
                            own,
                            nei
                        );

                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];
                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }
                        if (!modifiedFace)
                        {
                            modifiedFace = true;
                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }
            }

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }


    // 2. faces that do not get split but use edges that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Adding edge splits to unsplit faces"
            << endl;
    }

    DynamicList<label> eFacesStorage;
    DynamicList<label> fEdgesStorage;

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] >= 0)
        {
            // Split edge. Check that face not already handled above.

            const labelList& eFaces = mesh_.edgeFaces(edgeI, eFacesStorage);

            forAll(eFaces, i)
            {
                label facei = eFaces[i];

                if (faceMidPoint[facei] < 0 && affectedFace.get(facei))
                {
                    // Unsplit face. Add edge splits to face.

                    const face& f = mesh_.faces()[facei];
                    const labelList& fEdges = mesh_.faceEdges
                    (
                        facei,
                        fEdgesStorage
                    );

                    DynamicList<label> newFaceVerts(f.size());

                    forAll(f, fp)
                    {
                        newFaceVerts.append(f[fp]);

                        label edgeI = fEdges[fp];

                        if (edgeMidPoint[edgeI] >= 0)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    face newFace;
                    newFace.transfer(newFaceVerts);

                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    label anchorFp = findMinLevel(f);

                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        facei,
                        f[anchorFp],          // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        if (mesh_.isInternalFace(facei))
                        {
                            label oldOwn = mesh_.faceOwner()[facei];
                            label oldNei = mesh_.faceNeighbour()[facei];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                facei,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[facei];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                facei,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[facei],
                                newFace
                            );
                        }
                    }

                    modFace(meshMod, facei, newFace, own, nei);

                    // Mark face as having been handled
                    affectedFace.unset(facei);
                }
            }
        }
    }


    // 3. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }

    forAll(affectedFace, facei)
    {
        if (affectedFace.get(facei))
        {
            const face& f = mesh_.faces()[facei];

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            label anchorFp = findMinLevel(f);

            label own, nei;
            getFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                facei,
                f[anchorFp],          // Anchor point

                own,
                nei
            );

            modFace(meshMod, facei, f, own, nei);

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }


    // 4. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // This is the hard one. We have to find the splitting points between
    // the anchor points. But the edges between the anchor points might have
    // been split (into two,three or four edges).

    if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, celli)
    {
        if (cellMidPoint[celli] >= 0)
        {
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                celli,
                meshMod
            );
        }
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointi = labelMax;
        label maxPointi = labelMin;

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                minPointi = min(minPointi, cellMidPoint[celli]);
                maxPointi = max(maxPointi, cellMidPoint[celli]);
            }
        }
        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                minPointi = min(minPointi, faceMidPoint[facei]);
                maxPointi = max(maxPointi, faceMidPoint[facei]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointi = min(minPointi, edgeMidPoint[edgeI]);
                maxPointi = max(maxPointi, edgeMidPoint[edgeI]);
            }
        }

        if (minPointi != labelMax && minPointi != mesh_.nPoints())
        {
            FatalErrorInFunction
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointi:" << minPointi
                << " maxPointi:" << maxPointi
                << abort(FatalError);
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexRef4::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, celli)
        {
            const labelList& addedCells = cellAddedCells[celli];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(celli, addedCells);
            }
        }
    }

    // Compact cellAddedCells.

    labelListList refinedCells(cellLabels.size());

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        refinedCells[i].transfer(cellAddedCells[celli]);
    }

    return refinedCells;
}

Foam::labelList Foam::hexRef4::selectUnrefineElems
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitEdges(getSplitElems());

    DynamicList<label> newSplitEdges(splitEdges.size());

    forAll(splitEdges, i)
    {
        label edgej = splitEdges[i];

        const edge& e = mesh_.edges()[edgej];

        forAll(e, i)
        {
            label pointi = e[i];

            bool hasMarked = true;

            if (pFld[pointi] < unrefineLevel)
            {
                // Check that all cells are not marked
                const labelList& pCells = mesh_.pointCells()[pointi];

                hasMarked = false;

                forAll(pCells, pCelli)
                {
                    if (markedCell.get(pCells[pCelli]))
                    {
                        hasMarked = true;
                        break;
                    }
                }
            }

            if (!hasMarked)
            {
                newSplitEdges.append(edgej);
                break;
            }
        }
    }


    newSplitEdges.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        consistentUnrefinement
        (
            newSplitEdges,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split edges out of a possible "
        << returnReduce(splitEdges.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}

Foam::labelList Foam::hexRef4::consistentUnrefinement
(
    const labelList& elemsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef4::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorInFunction
            << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // For hexRef4, unrefinement is based on edges
    const labelList& edgesToUnrefine(elemsToUnrefine);

    // Loop, modifying edgesToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect edges to refine
    // maxSet = true: select edges to refine

    // Maintain boolList for edgesToUnrefine and cellsToUnrefine
    PackedBoolList unrefineEdge(mesh_.nEdges());

    forAll(edgesToUnrefine, i)
    {
        label edgei = edgesToUnrefine[i];

        unrefineEdge.set(edgei);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const labelListList& edgeCells = mesh_.edgeCells();

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefineEdge, edgei)
        {
            if (unrefineEdge.get(edgei))
            {
                const labelList& pCells = edgeCells[edgei];

                forAll(pCells, j)
                {
                    unrefineCell.set(pCells[j]);
                }
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            if (ownLevel < (neiLevel-1))
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    // could also combine with unset:
                    // if (!unrefineCell.unset(own))
                    // {
                    //     FatalErrorInFunction
                    //         << "problem cell already unset"
                    //         << abort(FatalError);
                    // }
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if (neiLevel < (ownLevel-1))
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            if (ownLevel < (neiLevel[i]-1))
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if (neiLevel[i] < (ownLevel-1))
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any edge whose cell neighbour cannot be unrefined.
        forAll(unrefineEdge, edgei)
        {
            if (unrefineEdge.get(edgei))
            {
                const labelList& pCells = edgeCells[edgei];

                forAll(pCells, j)
                {
                    if (!unrefineCell.get(pCells[j]))
                    {
                        unrefineEdge.unset(edgei);
                        break;
                    }
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefineEdge, edgei)
    {
        if (unrefineEdge.get(edgei))
        {
            nSet++;
        }
    }

    labelList newEdgesToUnrefine(nSet);
    nSet = 0;

    forAll(unrefineEdge, edgei)
    {
        if (unrefineEdge.get(edgei))
        {
            newEdgesToUnrefine[nSet++] = edgei;
        }
    }

    return newEdgesToUnrefine;
}

void Foam::hexRef4::calcFaceToSplitPoint
(
    const labelList& splitElems,
    Map<label>& faceToSplitPoint
)
{
    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // For hexRef4, the split elems are edges
    const labelList& splitEdges(splitElems);

    faceToSplitPoint.resize(2*splitEdges.size());

    {
        forAll(splitEdges, i)
        {
            label edgei = splitEdges[i];

            const edge& e = mesh_.edges()[edgei];

            forAll(e, j)
            {
                label pointi = e[j];

                const labelList& pFaces = mesh_.pointFaces()[pointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], pointi);
                }
            }
        }
    }
}

Foam::labelList Foam::hexRef4::getSplitElems() const
{
    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexRef4::getSplitElems :"
            << " Calculating unrefineable mid elements" << endl;
    }


    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split edge
    // >= label of master cell
    labelList splitMaster(mesh_.nEdges(), -1);
    labelList splitMasterLevel(mesh_.nEdges(), 0);

    // Unmark all with not 4 cells (so internal edges)
    const labelListList& edgeCells = mesh_.edgeCells();

    forAll(edgeCells, edgei)
    {
        const labelList& eCells = edgeCells[edgei];

        if (eCells.size() != 4)
        {
            splitMaster[edgei] = -2;
        }
    }

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, celli)
    {
        const labelList& cEdges = mesh_.cellEdges(celli);

        if (visibleCells[celli] != -1 && history_.parentIndex(celli) >= 0)
        {
            label parentIndex = history_.parentIndex(celli);

            // Check same master.
            forAll(cEdges, i)
            {
                label edgei = cEdges[i];

                label masterCelli = splitMaster[edgei];

                if (masterCelli == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to celli). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[edgei] = parentIndex;
                    splitMasterLevel[edgei] = cellLevel_[celli] - 1;
                }
                else if (masterCelli == -2)
                {
                    // Already decided that edge is not splitEdge
                }
                else if
                (
                    (masterCelli != parentIndex)
                 || (splitMasterLevel[edgei] != cellLevel_[celli] - 1)
                )
                {
                    // Different masters so edge is on two refinement
                    // patterns
                    splitMaster[edgei] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cEdges, i)
            {
                label edgei = cEdges[i];

                splitMaster[edgei] = -2;
            }
        }
    }

    // Unmark boundary faces
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        const labelList& fEdges = mesh_.faceEdges()[facei];

        forAll(fEdges, i)
        {
            label edgei = fEdges[i];

            splitMaster[edgei] = -2;
        }
    }


    // Collect into labelList

    label nSplitEdges = 0;

    forAll(splitMaster, edgei)
    {
        if (splitMaster[edgei] >= 0)
        {
            nSplitEdges++;
        }
    }

    labelList splitEdges(nSplitEdges);
    nSplitEdges = 0;

    forAll(splitMaster, edgei)
    {
        if (splitMaster[edgei] >= 0)
        {
            splitEdges[nSplitEdges++] = edgei;
        }
    }

    return splitEdges;
}

void Foam::hexRef4::setUnrefinement
(
    const labelList& splitElemLabels,
    polyTopoChange& meshMod
)
{
    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // For hexRef4, unrefinement is based on edges
    const labelList& splitEdgeLabels(splitElemLabels);


    //YO- This debug info would require writing edgeSets, which we do not have
    //if (debug)
    //{
    //    Pout<< "hexRef::setUnrefinement :"
    //        << " Checking initial mesh just to make sure" << endl;

    //    checkMesh();

    //    forAll(cellLevel_, celli)
    //    {
    //        if (cellLevel_[celli] < 0)
    //        {
    //            FatalErrorInFunction
    //                << "Illegal cell level " << cellLevel_[celli]
    //                << " for cell " << celli
    //                << abort(FatalError);
    //        }
    //    }


    //    // Write to sets.
    //    pointSet pSet(mesh_, "splitPoints", splitPointLabels);
    //    pSet.write();

    //    cellSet cSet(mesh_, "splitPointCells", splitPointLabels.size());

    //    forAll(splitPointLabels, i)
    //    {
    //        const labelList& pCells = mesh_.pointCells(splitPointLabels[i]);

    //        forAll(pCells, j)
    //        {
    //            cSet.insert(pCells[j]);
    //        }
    //    }
    //    cSet.write();

    //    Pout<< "hexRef::setRefinement : Dumping " << pSet.size()
    //        << " points and "
    //        << cSet.size() << " cells for unrefinement to" << nl
    //        << "    pointSet " << pSet.objectPath() << nl
    //        << "    cellSet " << cSet.objectPath()
    //        << endl;
    //}


    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        labelHashSet splitFaces(4*splitEdgeLabels.size());

        forAll(splitEdgeLabels, i)
        {
            const labelList& eFaces = mesh_.edgeFaces()[splitEdgeLabels[i]];

            forAll(eFaces, j)
            {
                splitFaces.insert(eFaces[j]);
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorInFunction
                << "Initial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitEdgeLabels, i)
    {
        label edgei = splitEdgeLabels[i];

        // Get original cell label

        const labelList& eCells = mesh_.edgeCells()[edgei];

        // Check
        if (eCells.size() != 4)
        {
            FatalErrorInFunction
                << "splitEdge " << edgei
                << " should have " << 4
                << " cells using it. It has " << eCells
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            label masterCelli = min(eCells);

            forAll(eCells, j)
            {
                label celli = eCells[j];

                label region = cellRegion[celli];

                if (region == -1)
                {
                    FatalErrorInFunction
                        << "Ininitial set of split edges to unrefine does not"
                        << " seem to be consistent or not mid edges"
                        << " of refined cells" << nl
                        << "cell:" << celli << " on splitEdge " << edgei
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCelli != cellRegionMaster[region])
                {
                    FatalErrorInFunction
                        << "cell:" << celli << " on splitEdge:" << edgei
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the edgeCells:" << eCells
                        << abort(FatalError);
                }
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Remove the n cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitEdgeLabels, i)
    {
        label edgei = splitEdgeLabels[i];

        const labelList& eCells = mesh_.edgeCells()[edgei];

        label masterCelli = min(eCells);

        forAll(eCells, j)
        {
            cellLevel_[eCells[j]]--;
        }

        history_.combineCells(masterCelli, eCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// ************************************************************************* //
