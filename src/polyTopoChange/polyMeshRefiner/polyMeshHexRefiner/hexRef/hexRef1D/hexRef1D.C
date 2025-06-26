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

#include "hexRef1D.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "meshTools.H"
#include "blastMeshTools.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef1D, 0);
    addToRunTimeSelectionTable(hexRef, hexRef1D, mesh);
    addToRunTimeSelectionTable(hexRef, hexRef1D, levelsHist);
    addToRunTimeSelectionTable(hexRef, hexRef1D, levels);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::hexRef1D::getAnchorCell
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
            // if (index >= 4) //AB....
            // {
            //     if (index == 4)
            //     {
            //         index = 8;
            //     }
            //     index = 8 - index;
            // } //AB
            return cellAddedCells[celli][index % 2];
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
                // if (index >= 4) //AB....
                // {
                //     if (index == 4)
                //     {
                //         index = 8;
                //     }
                //     index = 8 - index;
                // } //AB
                return cellAddedCells[celli][index % 2];
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


void Foam::hexRef1D::createInternalFace
(
    const labelListList& cellAddedCells,
    const HashTable<label, labelPair, Hash<labelPair>>& pointCellAnchorCell,
    const boolList& isEmptyFace,
    const labelList& edgeMidPoint,
    const label celli,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[celli];
    const labelList& cPoints = mesh_.cellPoints()[celli];
    const labelList& cEdges = mesh_.cellEdges()[celli];

    // Find faces on either side of the cell that are not empty
    label masterFace = -1;
    forAll(cFaces, i)
    {
        const label facei = cFaces[i];
        if (!isEmptyFace[facei])
        {
            masterFace = facei;
            break;
        }
    }

    const face& mf = mesh_.faces()[masterFace];

    // Storage for on-the-fly addressing
    DynamicList<label> storage(mf.size());

    forAll(mf, pi)
    {
        const label pointi = mf[pi];
        const labelList pEdges = mesh_.pointEdges()[pointi];
        forAll(pEdges, ei)
        {
            const label edgei = pEdges[ei];
            if (edgeMidPoint[edgei] >= 0 && findIndex(cEdges, edgei) >= 0)
            {
                storage.append(edgeMidPoint[edgei]);
                break;
            }
        }
    }

    label otherPoint = -1;
    forAll(cPoints, pi)
    {
        if (findIndex(mf, cPoints[pi]) < 0)
        {
            otherPoint = cPoints[pi];
            break;
        }
    }

    if (storage.size() != mf.size())
    {
        FatalErrorInFunction
            << "Could not split cell " << celli << endl
            << abort(FatalError);
    }

    label own = pointCellAnchorCell[{mf[0], celli}];
    label nei = pointCellAnchorCell[{otherPoint, celli}];

    face newFace(storage);
    if ((own != celli) == (mesh_.faceOwner()[masterFace] != celli))
    {
        newFace.flip();
    }

    meshTools::addInternalFace
    (
        meshMod,
        mesh_,
        masterFace,
        mf[0],
        newFace,
        own,
        nei

    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef1D::hexRef1D(const polyMesh& mesh, const bool readHistory)
:
    hexRef(mesh, readHistory)
{}


// Construct from components
Foam::hexRef1D::hexRef1D
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const hexRefRefinementHistory& history,
    const scalar level0Edge
)
:
    hexRef(mesh, cellLevel, pointLevel, history, level0Edge)
{}


// Construct from components
Foam::hexRef1D::hexRef1D
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
Foam::labelListList Foam::hexRef1D::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef1D::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
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

    locationMapper& locMapper(locationMapper::NewRef(mesh_));
    locMapper.clearOut();

    if (debug)
    {
        Pout<< "hexRef1D::setRefinement :"
            << " Allocating " << cellLabels.size() << " cell midpoints."
            << endl;
    }


    static const scalar angleTol =
        cos(Foam::constant::mathematical::pi/3.0);

    boolList isEmptyFace(mesh_.nFaces(), false);
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        if
        (
            isA<emptyPolyPatch>(patch)
         || isA<wedgePolyPatch>(patch)
        )
        {
            forAll(patch, fi)
            {
                isEmptyFace[patch.start() + fi] = true;
            }
        }
    }

    // Only mark edges that have 2 empty faces connected and an
    // angle greater than 45deg (should be ~90deg)
    boolList isDivisibleEdge(mesh_.nEdges(), false);
    const labelListList& edgeFaces = mesh_.edgeFaces();
    const scalarField& magFaceAreas = mesh_.magFaceAreas();
    const vectorField& faceAreas = mesh_.faceAreas();
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        if
        (
            isA<emptyPolyPatch>(patch)
         || isA<wedgePolyPatch>(patch)
        )
        {
            const labelList& meshEdges = patch.meshEdges();
            forAll(meshEdges, ei)
            {
                const label edgei = meshEdges[ei];
                const labelList& eFaces = edgeFaces[edgei];
                label face0 = -1;
                label face1 = -1;
                if (edgeFaces.size() < 2)
                {}
                else if (edgeFaces.size() == 2)
                {
                    face0 = isEmptyFace[eFaces[0]] ? eFaces[0] : -1;
                    face1 = isEmptyFace[eFaces[1]] ? eFaces[1] : -1;
                }
                else
                {
                    forAll(eFaces, fi)
                    {
                        const label facei = eFaces[fi];
                        if (isEmptyFace[facei])
                        {
                            if (face0 < 0)
                            {
                                face0 = facei;
                            }
                            else if (face1 < 0)
                            {
                                face1 = facei;
                            }
                            else
                            {
                                FatalErrorInFunction
                                    << "Edge connected to more than 2 "
                                    << "empty faces" << endl;
                            }
                        }
                    }
                }

                // Check angle between faces
                if (face0 >= 0 && face1 >= 0)
                {
                    const vector nf0 =
                        faceAreas[face0]/magFaceAreas[face0];
                    const vector nf1 =
                        faceAreas[face1]/magFaceAreas[face1];

                    if (mag(nf0 & nf1) < angleTol)
                    {
                        isDivisibleEdge[edgei] = true;
                    }
                }
            }
        }
    }


    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellLabels);
        Pout<< "hexRef1D::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }



    // Split edges
    // ~~~~~~~~~~~
    DebugInFunction<< "Allocating edge midpoints" << endl;

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(edges.size(), -1);

    // Note: Loop over cells to be refined or edges?
    forAll(cellLabels, ci)
    {
        const label celli = cellLabels[ci];
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
                edgeMidPoint[edgeI] = 12345;
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

        // Add split edges
        labelList newEdgePoints(edgeMidPoint.size(), -1);
        pointField edgeMids(edges.size(), point::zero);

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = edges[edgeI].centre(points);
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxMagSqrEqOp<vector>(),
            point::zero
        );


        // Phase 2: introduce points at the synced locations.
        DynamicList<label> splitEdges(edgeMidPoint.size());
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.
                const edge& e = edges[edgeI];

                edgeMidPoint[edgeI] = meshMod.addPoint
                (
                    edgeMids[edgeI],            // point
                    e[0],                       // master point
                    true                        // supports a cell
                );
                splitEdges.append(edgeI);
                newEdgePoints[edgeI] = edgeMidPoint[edgeI];

                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    ) + 1;
            }
        }
        locMapper.addSplitEdges(splitEdges, newEdgePoints);

        if (debug)
        {
            OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

            forAll(edgeMidPoint, edgeI)
            {
                if (edgeMidPoint[edgeI] >= 0)
                {
                    meshTools::writeOBJ(str, edgeMids[edgeI]);
                }
            }

            Pout<< FUNCTION_NAME << ": "
                << "Dumping edge centres to split to file "
                << str.name() << endl;
        }
    }



    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - edgeMidPoint >= 0 : edge needs to be split.
    // - isEmptyFace true : face needs to be split in 2


    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef1D::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());
    {
        labelList nAnchorPoints(mesh_.nCells(), 0);
        forAll(cellLabels, ci)
        {
            cellAnchorPoints[cellLabels[ci]].setSize(8);
        }

       forAll(cellLabels, ci)
       {
            const label celli = cellLabels[ci];
            const labelList& cPoints = mesh_.cellPoints()[celli];
            forAll(cPoints, i)
            {
                const label pointi = cPoints[i];
                if (pointLevel_[pointi] <= cellLevel_[celli])
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


        forAll(cellLabels, ci)
        {
            const label celli = cellLabels[ci];
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


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< FUNCTION_NAME << ": "
            << " Adding cells (1 per anchorPoint)"
            << endl;
    }

    // Per cell the 1 added cells (+ original cell)
    labelListList cellAddedCells(mesh_.nCells());
    HashTable<label, labelPair, Hash<labelPair>> pointCellAnchorCell;

    forAll(cellAnchorPoints, celli)
    {
        const labelList& cAnchors = cellAnchorPoints[celli];

        if (cAnchors.size() == 8)
        {
            labelList& cAdded = cellAddedCells[celli];
            cAdded.setSize(2);

            // Original cell at 0
            cAdded[0] = celli;
            cAdded[1] = meshMod.addCell(celli);

            // Update cell levels
            newCellLevel[celli] = cellLevel_[celli]+1;
            newCellLevel(cAdded[1]) = cellLevel_[celli]+1;

            label masterFace = -1;
            const cell& c = mesh_.cells()[celli];
            forAll(c, fi)
            {
                if (!isEmptyFace[c[fi]])
                {
                    masterFace = c[fi];
                    break;
                }
            }

            const face& mf = mesh_.faces()[masterFace];
            const labelList& cPoints = mesh_.cellPoints()[celli];
            forAll(cPoints, pi)
            {
                const label pointi = cPoints[pi];
                if (findIndex(mf, pointi) >= 0)
                {
                    pointCellAnchorCell.insert
                    (
                        labelPair(pointi, celli),
                        cAdded[0]
                    );
                }
                else
                {
                    pointCellAnchorCell.insert
                    (
                        labelPair(pointi, celli),
                        cAdded[1]
                    );
                }
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into two always)
    // 3. existing faces that do not get split but get new owner/neighbour
    // 4. new internal faces inside split cells.

    DebugInFunction<< " Marking faces to be handled" << endl;

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    forAll(cellLabels, ci)
    {
        const cell& cFaces = mesh_.cells()[cellLabels[ci]];

        forAll(cFaces, i)
        {
            affectedFace.set(cFaces[i]);
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef1D::setRefinement : Splitting faces" << endl;
    }

    forAll(affectedFace, facei)
    {
        if (affectedFace.get(facei) && isEmptyFace[facei])
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const label celli = mesh_.faceOwner()[facei];

            const face& f = mesh_.faces()[facei];

            // Has original facei been used (three faces added, original gets
            // modified)
            const label anchorLevel = faceLevel(facei);

            // Find start point of walk
            label pi = 0;
            forAll(f, i)
            {
                if (pointLevel_[f[i]] <= anchorLevel)
                {
                    pi = i;
                    break;
                }
            }

            bool isFace0 = true;
            DynamicList<label> faceVerts0(4), faceVerts1(4);
            forAll(f, i)
            {
                const label pointi = f[pi];
                const label pointj = f[f.fcIndex(pi)];
                const label edgei = meshTools::findEdge(mesh_, pointi, pointj);

                // point is anchor. Start collecting face.

                if (isFace0)
                {
                    faceVerts0.append(pointi);
                }
                else
                {
                    faceVerts1.append(pointi);
                }

                if (edgeMidPoint[edgei] >= 0)
                {
                    faceVerts0.append(edgeMidPoint[edgei]);
                    faceVerts1.append(edgeMidPoint[edgei]);
                    isFace0 = !isFace0;
                }

                pi++;
            }

            face newFace0(faceVerts0);
            face newFace1(faceVerts1);
            if (debug)
            {
                meshTools::checkFaceOrientation
                (
                    meshMod,
                    mesh_,
                    facei,
                    newFace0
                );
                meshTools::checkFaceOrientation
                (
                    meshMod,
                    mesh_,
                    facei,
                    newFace1
                );
            }
            {
                meshTools::modifyFace
                (
                    meshMod,
                    mesh_,
                    facei,
                    newFace0,
                    pointCellAnchorCell[{min(newFace0), celli}],
                    -1
                );
                meshTools::addFace
                (
                    meshMod,
                    mesh_,
                    facei,
                    newFace1,
                    pointCellAnchorCell[{min(newFace1), celli}],
                    -1
                );
            }

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }

    // 2. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DebugInFunction
        << " Changing owner/neighbour for otherwise unaffected faces"
        << endl;

    forAll(affectedFace, facei)
    {
        if (affectedFace.get(facei))
        {
            const face& f = mesh_.faces()[facei];

            label own = mesh_.faceOwner()[facei];
            label nei =
                facei < mesh_.nInternalFaces()
              ? mesh_.faceNeighbour()[facei]
              : -1;

            face newFace(f);

            // Correct own/nei to use new cell indices
            HashTable
            <
                label,
                labelPair,
                Hash<labelPair>
            >::const_iterator iterOwn = pointCellAnchorCell.find({f[0], own});
            if (iterOwn != pointCellAnchorCell.cend())
            {
                own = iterOwn();
            }

            if (nei >= 0)
            {
                HashTable
                <
                    label,
                    labelPair,
                    Hash<labelPair>
                >::const_iterator iterNei =
                    pointCellAnchorCell.find({f[0], nei});
                if (iterNei != pointCellAnchorCell.cend())
                {
                    nei = iterNei();
                }

                if (own > nei)
                {
                    Swap(own, nei);
                    newFace.flip();
                }
            }

            // Modify exisiting face
            meshTools::modifyFace
            (
                meshMod,
                mesh_,
                facei,
                newFace,
                own,
                nei
            );

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
        Pout<< "hexRef1D::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellLabels, ci)
    {
        createInternalFace
        (
            cellAddedCells,
            pointCellAnchorCell,
            isEmptyFace,
            edgeMidPoint,
            cellLabels[ci],
            meshMod
        );
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointi = labelMax;
        label maxPointi = labelMin;
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
            Pout<< "hexRef1D::setRefinement :"
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

Foam::labelList Foam::hexRef1D::selectUnrefineElems
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitFaces(getSplitElems());

    DynamicList<label> newSplitFaces(splitFaces.size());

    forAll(splitFaces, i)
    {
        label facej = splitFaces[i];

        const face& f = mesh_.faces()[facej];

        forAll(f, pi)
        {
            label pointi = f[pi];

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
                newSplitFaces.append(facej);
                break;
            }
        }
    }
    newSplitFaces.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        consistentUnrefinement
        (
            newSplitFaces,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split faces out of a possible "
        << returnReduce(splitFaces.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}

Foam::labelList Foam::hexRef1D::consistentUnrefinement
(
    const labelList& elemsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef1D::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorInFunction
            << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // For hexRef1D, unrefinement is based on edges
    const labelList& facesToUnrefine(elemsToUnrefine);

    // Loop, modifying edgesToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect edges to refine
    // maxSet = true: select edges to refine

    // Maintain boolList for edgesToUnrefine and cellsToUnrefine
    PackedBoolList unrefineFace(mesh_.nInternalFaces());

    forAll(facesToUnrefine, i)
    {
        label facei = facesToUnrefine[i];

        unrefineFace.set(facei);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefineFace, facei)
        {
            if (unrefineFace.get(facei))
            {
                const label own = mesh_.faceOwner()[facei];
                const label nei = mesh_.faceNeighbour()[facei];

                unrefineCell.set(own);
                unrefineCell.set(nei);
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

        DebugInFunction
            << " Changed " << nChanged
            << " refinement levels due to 2:1 conflicts."
            << endl;

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any face whose cell neighbour cannot be unrefined.
        forAll(unrefineFace, facei)
        {
            if (unrefineFace.get(facei))
            {
                const label own = mesh_.faceOwner()[facei];
                const label nei = mesh_.faceNeighbour()[facei];
                if (!unrefineCell.get(own) || !unrefineCell.get(nei))
                {
                    unrefineFace.unset(facei);
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefineFace, facei)
    {
        if (unrefineFace.get(facei))
        {
            nSet++;
        }
    }

    labelList newFacesToUnrefine(nSet);
    nSet = 0;

    forAll(unrefineFace, facei)
    {
        if (unrefineFace.get(facei))
        {
            newFacesToUnrefine[nSet++] = facei;
        }
    }

    return newFacesToUnrefine;
}

void Foam::hexRef1D::calcFaceToSplitPoint
(
    const labelList& splitElems,
    Map<label>& faceToSplitPoint
)
{
    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // For hexRef1D, the split elems are faces
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

Foam::labelList Foam::hexRef1D::getSplitElems() const
{
    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexRef1D::getSplitElems :"
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
    // -2 certainly not split face
    // >= label of master cell
    labelList splitMaster(mesh_.nInternalFaces(), -1);
    labelList splitMasterLevel(mesh_.nInternalFaces(), 0);

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, celli)
    {
        const labelList& c = mesh_.cells()[celli];

        if (visibleCells[celli] != -1 && history_.parentIndex(celli) >= 0)
        {
            label parentIndex = history_.parentIndex(celli);

            // Check same master.
            forAll(c, fi)
            {
                const label facei = c[fi];

                if (!mesh_.isInternalFace(facei))
                {
                    continue;
                }

                label masterCelli = splitMaster[facei];

                if (masterCelli == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to celli). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[facei] = parentIndex;
                    splitMasterLevel[facei] = cellLevel_[celli] - 1;
                }
                else if (masterCelli == -2)
                {
                    // Already decided that face is not splitFace
                }
                else if
                (
                    (masterCelli != parentIndex)
                 || (splitMasterLevel[facei] != cellLevel_[celli] - 1)
                )
                {
                    // Different masters so face is on two refinement
                    // patterns
                    splitMaster[facei] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(c, fi)
            {
                if (mesh_.isInternalFace(c[fi]))
                {
                    splitMaster[c[fi]] = -2;
                }
            }
        }
    }


    // Collect into labelList

    // Count split faces
    label nSplitFaces = 0;
    forAll(splitMaster, facei)
    {
        if (splitMaster[facei] >= 0)
        {
            nSplitFaces++;
        }
    }

    // Insert split faces
    labelList splitFaces(nSplitFaces);
    nSplitFaces = 0;
    forAll(splitMaster, facei)
    {
        if (splitMaster[facei] >= 0)
        {
            splitFaces[nSplitFaces++] = facei;
        }
    }

    return splitFaces;
}

void Foam::hexRef1D::setUnrefinement
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

    // For hexRef1D, unrefinement is based on faces
    const labelList& splitFaceLabels(splitElemLabels);


    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;
    {
        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaceLabels,    // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaceLabels.size())
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

    forAll(splitFaceLabels, i)
    {
        const label facei = splitFaceLabels[i];

        // Check
        if (!mesh_.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "splitFace " << facei
                << " should be internal but it is not" << endl
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        if (debug)
        {
            const label own = mesh_.faceOwner()[facei];
            const label nei = mesh_.faceNeighbour()[facei];
            label masterCelli = min(own, nei);

            {
                const label ownRegion = cellRegion[own];
                const label neiRegion = cellRegion[own];

                if (ownRegion == -1)
                {
                    FatalErrorInFunction
                        << "Ininitial set of split faces to unrefine does not"
                        << " seem to be consistent" << nl
                        << "cell:" << own << " on splitFace " << facei
                        << " has no region to be merged into"
                        << abort(FatalError);
                }
                if (neiRegion == -1)
                {
                    FatalErrorInFunction
                        << "Ininitial set of split faces to unrefine does not"
                        << " seem to be consistent" << nl
                        << "cell:" << nei << " on splitFace " << facei
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCelli != cellRegionMaster[ownRegion])
                {
                    FatalErrorInFunction
                        << "cell:" << own << " on splitFace:" << facei
                        << " in region " << ownRegion
                        << " has master:" << cellRegionMaster[ownRegion]
                        << " which is not the lowest numbered cell"
                        << " among the own/nei:" << own << "/" << nei << endl
                        << abort(FatalError);
                }
                if (masterCelli != cellRegionMaster[neiRegion])
                {
                    FatalErrorInFunction
                        << "cell:" << nei << " on splitFace:" << facei
                        << " in region " << neiRegion
                        << " has master:" << cellRegionMaster[neiRegion]
                        << " which is not the lowest numbered cell"
                        << " among the own/nei:" << own << "/" << nei << endl
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
    labelList fCells(2);
    forAll(splitFaceLabels, i)
    {
        const label facei = splitFaceLabels[i];

        fCells = {mesh_.faceOwner()[facei], mesh_.faceNeighbour()[facei]};

        label masterCelli = min(fCells);

        forAll(fCells, j)
        {
            cellLevel_[fCells[j]]--;
        }

        history_.combineCells(masterCelli, fCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// ************************************************************************* //
