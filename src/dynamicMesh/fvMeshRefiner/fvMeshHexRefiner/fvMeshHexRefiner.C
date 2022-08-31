/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
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

#include "fvMeshHexRefiner.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"
#include "cellSet.H"
#include "wedgePolyPatch.H"
#include "hexRef3D.H"
#include "parcelCloud.H"
#include "hexRefRefinementHistoryConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshHexRefiner, 0);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshHexRefiner, fvMesh);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshHexRefiner, dictionary);

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshHexRefiner::calculateProtectedCells
(
    PackedBoolList& unrefineableCell
) const
{
    if (!returnReduce(protectedCell_.size(), sumOp<label>()))
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_->cellLevel();

    unrefineableCell = protectedCell_;

    // Get neighbouring cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        neiLevel[facei-mesh_.nInternalFaces()] =
            cellLevel[mesh_.faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(mesh_.nFaces(), false);

        forAll(mesh_.faceNeighbour(), facei)
        {
            label own = mesh_.faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            label nei = mesh_.faceNeighbour()[facei];
            bool neiProtected = unrefineableCell.get(nei);

            if (ownProtected && (cellLevel[nei] > cellLevel[own]))
            {
                seedFace[facei] = true;
            }
            else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
            {
                seedFace[facei] = true;
            }
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            label own = mesh_.faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            if
            (
                ownProtected
             && (neiLevel[facei-mesh_.nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh_, seedFace, orEqOp<bool>());


        // Extend unrefineableCell
        bool hasExtended = false;

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = mesh_.faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }

                label nei = mesh_.faceNeighbour()[facei];
                if (unrefineableCell.get(nei) == 0)
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            if (seedFace[facei])
            {
                label own = mesh_.faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }
            }
        }

        if (!returnReduce(hasExtended, orOp<bool>()))
        {
            break;
        }
    }
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::fvMeshHexRefiner::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh_);

    // Play refinement commands into mesh changer.
    meshCutter_->setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh_.globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label oldFacei = map().faceMap()[facei];

            if (oldFacei >= mesh_.nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << mesh_.faceCentres()[facei]
                    << " originates from boundary oldFace:" << oldFacei
                    << abort(FatalError);
            }
        }
    }

    //    // Remove the stored tet base points
    //    tetBasePtIsPtr_.clear();
    //    // Remove the cell tree
    //    cellTreePtr_.clear();

    // Update fields
    mesh_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(mesh_.nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            newProtectedCell.set(celli, protectedCell_.get(oldCelli));
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_->checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh>
Foam::fvMeshHexRefiner::unrefine
(
    const labelList& splitElems
)
{
    polyTopoChange meshMod(mesh_);

    // Play refinement commands into mesh changer.
    meshCutter_->setUnrefinement(splitElems, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(0);
    meshCutter_->calcFaceToSplitPoint(splitElems, faceToSplitPoint);


    // Change mesh and generate map.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh_.globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    mesh_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(mesh_.nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            if (oldCelli >= 0)
            {
                newProtectedCell.set(celli, protectedCell_.get(oldCelli));
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_->checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::labelList Foam::fvMeshHexRefiner::selectRefineCells
(
    const label maxCells,
    const labelList& maxRefinement,
    const PackedBoolList& candidateCell
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - mesh_.globalData().nTotalCells()) / 7;
    const labelList& cellLevel = meshCutter_->cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCell;
    calculateProtectedCells(unrefineableCell);

    // Count current selection
    label nLocalCandidates = count(candidateCell, 1);
    label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);
    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCell, celli)
        {
            if
            (
                cellLevel[celli] < maxRefinement[celli]
             && candidateCell.get(celli)
             && (
                    unrefineableCell.empty()
                 || !unrefineableCell.get(celli)
                )
            )
            {
                candidates.append(celli);
            }
        }
    }
    else
    {
        label maxLevel(gMax(maxRefinement));

        // Sort by error? For now just truncate.
        for (label level = 0; level < maxLevel; level++)
        {
            forAll(candidateCell, celli)
            {
                if
                (
                    cellLevel[celli] == level
                 && candidateCell.get(celli)
                 && (
                        unrefineableCell.empty()
                     || !unrefineableCell.get(celli)
                    )
                )
                {
                    candidates.append(celli);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_->consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << mesh_.globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


void Foam::fvMeshHexRefiner::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    labelList nAnchorPoints(mesh_.nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = mesh_.pointCells(pointi);

        forAll(pCells, pCelli)
        {
            label celli = pCells[pCelli];

            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[celli] == 8)
                {
                    if (protectedCell.set(celli, true))
                    {
                        nProtected++;
                    }
                }

                if (!protectedCell.get(celli))
                {
                    nAnchorPoints[celli]++;
                }
            }
        }
    }


    forAll(protectedCell, celli)
    {
        if (!protectedCell.get(celli) && nAnchorPoints[celli] != 8)
        {
            protectedCell.set(celli, true);
            nProtected++;
        }
    }
}


void Foam::fvMeshHexRefiner::setProtectedCells()
{
    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    labelList neiLevel(mesh_.nFaces());

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        neiLevel[facei] = cellLevel[mesh_.faceNeighbour()[facei]];
    }
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        neiLevel[facei] = cellLevel[mesh_.faceOwner()[facei]];
    }
    syncTools::swapFaceList(mesh_, neiLevel);

    // Protect cells that will cause a failure (from snappyHexMesh)
    boolList protectedFaces(mesh_.nFaces(), false);

    forAll(mesh_.faceOwner(), facei)
    {
        label faceLevel = max
        (
            cellLevel[mesh_.faceOwner()[facei]],
            neiLevel[facei]
        );

        const face& f = mesh_.faces()[facei];

        label nAnchors = 0;

        forAll(f, fp)
        {
            if (pointLevel[f[fp]] <= faceLevel)
            {
                nAnchors++;
            }
        }
        if (nAnchors == 2)
        {
            protectedFaces[facei] = true;
        }
    }

    syncTools::syncFaceList(mesh_, protectedFaces, orEqOp<bool>());

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (protectedFaces[facei])
        {
            protectedCell_.set(mesh_.faceOwner()[facei], 1);
            nProtected_++;
            protectedCell_.set(mesh_.faceNeighbour()[facei], 1);
            nProtected_++;
        }
    }
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        if (protectedFaces[facei])
        {
            protectedCell_.set(mesh_.faceOwner()[facei], 1);
            nProtected_++;
        }
    }

    if (mesh_.nGeometricD() == mesh_.nSolutionD())
    {
        checkEightAnchorPoints(protectedCell_, nProtected_);
    }

    if (returnReduce(nProtected_, sumOp<label>()) == 0)
    {
        protectedCell_.clear();
    }
}


void Foam::fvMeshHexRefiner::updateMesh(const mapPolyMesh& mpm)
{
    fvMeshRefiner::updateMesh(mpm);

    // Do not update hexMesh since it is handled in the distribute function
    if (!isBalancing_)
    {
        meshCutter_->updateMesh(mpm);
    }
}


void Foam::fvMeshHexRefiner::distribute
(
    const mapDistributePolyMesh& map
)
{
    fvMeshRefiner::distribute(map);

    meshCutter_->distribute(map);
    if (returnReduce(nProtected_, sumOp<label>()) > 0)
    {
        boolList protectedCell(protectedCell_.size(), false);
        forAll(protectedCell, i)
        {
            if (protectedCell_.get(i))
            {
                protectedCell[i] = true;
            }
        }
        map.distributeCellData(protectedCell);
        protectedCell_ = protectedCell;
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshHexRefiner::fvMeshHexRefiner(fvMesh& mesh)
:
    fvMeshRefiner(mesh),

    meshCutter_(hexRef::New(mesh_)),

    nProtected_(0),
    protectedCell_(mesh_.nCells(), 0)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            hexRefRefinementHistoryConstraint::typeName
        );
        balancer_.addConstraint("refinementHistory", refinementHistoryDict);
    }

    nProtected_ = 0;

    forAll(protectedCell_, celli)
    {
        if (protectedCell_.get(celli))
        {
            nProtected_++;
        }
    }


    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(mesh_.nCells(), 0);

    forAll(mesh_.pointCells(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];

        forAll(pCells, i)
        {
            label celli = pCells[i];

            if (!protectedCell_.get(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    nAnchors[celli]++;

                    if (nAnchors[celli] > 8)
                    {
                        protectedCell_.set(celli, 1);
                        nProtected_++;
                    }
                }
            }
        }
    }

    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(mesh_.nFaces());

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[mesh_.faceNeighbour()[facei]];
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            neiLevel[facei] = cellLevel[mesh_.faceOwner()[facei]];
        }
        syncTools::swapFaceList(mesh_, neiLevel);


        boolList protectedFace(mesh_.nFaces(), false);

        forAll(mesh_.faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[mesh_.faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = mesh_.faces()[facei];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[facei] = true;
                        break;
                    }
                }
            }

        }

        syncTools::syncFaceList(mesh_, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(mesh_.faceOwner()[facei], 1);
                nProtected_++;
                protectedCell_.set(mesh_.faceNeighbour()[facei], 1);
                nProtected_++;
            }
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(mesh_.faceOwner()[facei], 1);
                nProtected_++;
            }
        }

        bool isAxisym = (mesh_.nGeometricD() == 2 && mesh_.nSolutionD() == 3);

        if (!isAxisym)
        {
            // Also protect any cells that are less than hex
            forAll(mesh_.cells(), celli)
            {
                const cell& cFaces = mesh_.cells()[celli];

                if (cFaces.size() < 6)
                {
                    if (protectedCell_.set(celli, 1))
                    {
                        nProtected_++;
                    }
                }
                else
                {
                    forAll(cFaces, cFacei)
                    {
                        if (mesh_.faces()[cFaces[cFacei]].size() < 4)
                        {
                            if (protectedCell_.set(celli, 1))
                            {
                                nProtected_++;
                            }
                            break;
                        }
                    }
                }
            }

            // Check cells for 8 corner points
            checkEightAnchorPoints(protectedCell_, nProtected_);
        }
        else
        {
            PackedBoolList cellIsAxisPrism(mesh_.nCells(), false);
            label nAxisPrims = 0;

            // Do not protect prisms on the axis
            forAll(mesh_.cells(), celli)
            {
                const cell& cFaces = mesh_.cells()[celli];

                if (cFaces.size() == 5)
                {
                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];
                        if
                        (
                            !mesh_.isInternalFace(facei)
                         && mesh_.faces()[facei].size() != 4
                        )
                        {
                            label patchi =
                                mesh_.boundaryMesh().whichPatch(facei);
                            if
                            (
                                isA<wedgePolyPatch>
                                (
                                    mesh_.boundaryMesh()[patchi]
                                )
                            )
                            {
                                if (protectedCell_.set(celli, 1))
                                {
                                    nProtected_++;
                                    break;
                                }
                            }
                        }
                    }
                    if (!protectedCell_.get(celli))
                    {
                        cellIsAxisPrism.set(celli, 1);
                        nAxisPrims++;
                    }
                }
                else if (cFaces.size() < 5)
                {
                    if (protectedCell_.set(celli, 1))
                    {
                        nProtected_++;
                    }
                }
            }

            // Check cells for 8 corner points
            checkEightAnchorPoints(protectedCell_, nProtected_);

            // Unprotect prism cells on the axis
            protectedCell_ -= cellIsAxisPrism;
            nProtected_ -= nAxisPrims;



            // YO
            // Unprotect cells that have a refinement history
            const labelList& visibleCells =
                meshCutter_->history().visibleCells();
            forAll(visibleCells, celli)
            {
                if (visibleCells[celli] >= 0)
                {
                    if (protectedCell_.get(celli))
                    {
                        protectedCell_.unset(celli);
                        nProtected_--;
                    }
                }
            }

            // Look for former prisms cells on level 0. Since they do not have
            // a refinement history, they will be incorrectly protected.
            // In order to find former prisms, look for protectedCells with
            // exactly 2 points on the axis. If the number of refined
            // neighbours matches with the number of faces in excess of 5,
            // the cell is (most probably) a former prism.

            {
                // First, look for a wedge polyPatch, in order to get the wedge
                // plane normal and axis.
                vector axis(0,0,0);
                vector centreNormal(0,0,0);
                bool foundWedge = false;

                forAll(mesh_.boundaryMesh(), patchi)
                {
                    if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchi]))
                    {
                        foundWedge = true;

                        const wedgePolyPatch& wedgePatch =
                            refCast<const wedgePolyPatch>
                            (
                                mesh_.boundaryMesh()[patchi]
                            );

                        axis = wedgePatch.axis();
                        centreNormal = wedgePatch.centreNormal();
                        break;
                    }
                }

                if (!foundWedge)
                {
                    FatalErrorInFunction
                        << "Number of geometric dimensions and solution dimensions "
                        << "correspond to an axisymmetric case, but no wedgePolyPatch "
                        << " found in boundaryMesh." << nl
                        << exit(FatalError);
                }

                // Get radial direction on the mesh
                vector radialVector = axis ^ centreNormal;

                // Translate as a component. Since wedge meshes should be constrained
                // around one of the planes XY, XZ or YZ, we can directly get the
                // distance to the axis

                direction dir = 0;
                scalar maxComp = mag(radialVector[0]);

                if (mag(radialVector[1]) > maxComp)
                {
                    dir = 1;
                    maxComp = mag(radialVector[1]);
                }
                if (mag(radialVector[2]) > maxComp)
                {
                    dir = 2;
                }

                // Use the level0Edge length to calculate a safety margin for
                // detecting points on the axis.
                // We should be able to rely on level0, since we only check
                // unrefined cells.
                scalar tolerance = meshCutter_->level0EdgeLength() * 1e-6;

                const labelListList& cellPoints = mesh_.cellPoints();

                // Get max level across faces
                labelList maxFaceLevel(mesh_.nFaces());
                const labelList& cellLevel = meshCutter_->cellLevel();

                forAll(maxFaceLevel, facei)
                {
                    maxFaceLevel[facei] = cellLevel[mesh_.faceOwner()[facei]];
                }

                forAll(mesh_.faceNeighbour(), facei)
                {
                    maxFaceLevel[facei] =
                        max
                        (
                            maxFaceLevel[facei],
                            cellLevel[mesh_.faceNeighbour()[facei]]
                        );
                }
                syncTools::syncFaceList
                (
                    mesh_,
                    maxFaceLevel,
                    combineMaxOp<label>()
                );

                forAll(protectedCell_, celli)
                {
                    if (protectedCell_.get(celli))
                    {
                        // Check if the cell has exactly 2 points on the axis
                        label numPointsOnAxis = 0;
                        const labelList& cPoints = cellPoints[celli];
                        forAll(cPoints, pointi)
                        {
                            if
                            (
                                mag(mesh_.points()[cPoints[pointi]][dir])
                              < tolerance
                            )
                            {
                                numPointsOnAxis++;
                            }
                        }

                        if (numPointsOnAxis == 2)
                        {
                            // Calculate number of higher level neighbour cells.
                            // Use cellFaces instead of cellCells to deal with
                            // processor patches.
                            label numNeighboursWithHigherLevel = 0;

                            const labelList& cellFaces = mesh_.cells()[celli];
                            forAll(cellFaces, cellFacei)
                            {
                                if
                                (
                                    cellLevel[celli]
                                  < maxFaceLevel[cellFaces[cellFacei]]
                                )
                                {
                                    numNeighboursWithHigherLevel++;
                                }
                            }

                            // If the cell was a neighbour to refined cells, then each
                            // pair of refined neighbours introduces an additional face
                            if
                            (
                                2*(cellFaces.size() - 5)
                             == numNeighboursWithHigherLevel
                            )
                            {
                                protectedCell_.unset(celli);
                                nProtected_--;
                            }
                        }
                    }
                }
            }
        }
    }
    setProtectedCells();

    if (returnReduce(nProtected_, sumOp<label>()))
    {
        cellSet protectedCells(mesh_, "protectedCells", nProtected_);
        forAll(protectedCell_, celli)
        {
            if (protectedCell_.get(celli))
            {
                protectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected_, sumOp<label>())
            << " cells that are protected from refinement." << endl;
    }
}


Foam::fvMeshHexRefiner::fvMeshHexRefiner
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
:
    fvMeshRefiner(mesh, dict, force, read),

    meshCutter_
    (
        dict.lookupOrDefault<bool>("force3D", false)
      ? new hexRef3D(mesh_, read)
      : hexRef::New(mesh_, read).ptr()
    ),

    nProtected_(0),
    protectedCell_(mesh_.nCells(), 0)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict;
        refinementHistoryDict.add
        (
            "type",
            hexRefRefinementHistoryConstraint::typeName
        );
        balancer_.addConstraint("refinementHistory", refinementHistoryDict);
    }

    // Read static part of dictionary
    readDict(dict);

    nProtected_ = 0;

    forAll(protectedCell_, celli)
    {
        if (protectedCell_.get(celli))
        {
            nProtected_++;
        }
    }


    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(mesh_.nCells(), 0);

    forAll(mesh_.pointCells(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];

        forAll(pCells, i)
        {
            label celli = pCells[i];

            if (!protectedCell_.get(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    nAnchors[celli]++;

                    if (nAnchors[celli] > 8)
                    {
                        protectedCell_.set(celli, 1);
                        nProtected_++;
                    }
                }
            }
        }
    }

    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(mesh_.nFaces());

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[mesh_.faceNeighbour()[facei]];
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            neiLevel[facei] = cellLevel[mesh_.faceOwner()[facei]];
        }
        syncTools::swapFaceList(mesh_, neiLevel);


        boolList protectedFace(mesh_.nFaces(), false);

        forAll(mesh_.faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[mesh_.faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = mesh_.faces()[facei];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[facei] = true;
                        break;
                    }
                }
            }

        }

        syncTools::syncFaceList(mesh_, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(mesh_.faceOwner()[facei], 1);
                nProtected_++;
                protectedCell_.set(mesh_.faceNeighbour()[facei], 1);
                nProtected_++;
            }
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(mesh_.faceOwner()[facei], 1);
                nProtected_++;
            }
        }

        bool isAxisym = (mesh_.nGeometricD() == 2 && mesh_.nSolutionD() == 3);

        if (!isAxisym)
        {
            // Also protect any cells that are less than hex
            forAll(mesh_.cells(), celli)
            {
                const cell& cFaces = mesh_.cells()[celli];

                if (cFaces.size() < 6)
                {
                    if (protectedCell_.set(celli, 1))
                    {
                        nProtected_++;
                    }
                }
                else
                {
                    forAll(cFaces, cFacei)
                    {
                        if (mesh_.faces()[cFaces[cFacei]].size() < 4)
                        {
                            if (protectedCell_.set(celli, 1))
                            {
                                nProtected_++;
                            }
                            break;
                        }
                    }
                }
            }

            // Check cells for 8 corner points
            checkEightAnchorPoints(protectedCell_, nProtected_);
        }
        else
        {
            PackedBoolList cellIsAxisPrism(mesh_.nCells(), false);
            label nAxisPrims = 0;

            // Do not protect prisms on the axis
            forAll(mesh_.cells(), celli)
            {
                const cell& cFaces = mesh_.cells()[celli];

                if (cFaces.size() == 5)
                {
                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];
                        if
                        (
                            !mesh_.isInternalFace(facei)
                         && mesh_.faces()[facei].size() != 4
                        )
                        {
                            label patchi =
                                mesh_.boundaryMesh().whichPatch(facei);
                            if
                            (
                                isA<wedgePolyPatch>
                                (
                                    mesh_.boundaryMesh()[patchi]
                                )
                            )
                            {
                                if (protectedCell_.set(celli, 1))
                                {
                                    nProtected_++;
                                    break;
                                }
                            }
                        }
                    }
                    if (!protectedCell_.get(celli))
                    {
                        cellIsAxisPrism.set(celli, 1);
                        nAxisPrims++;
                    }
                }
                else if (cFaces.size() < 5)
                {
                    if (protectedCell_.set(celli, 1))
                    {
                        nProtected_++;
                    }
                }
            }

            // Check cells for 8 corner points
            checkEightAnchorPoints(protectedCell_, nProtected_);

            // Unprotect prism cells on the axis
            protectedCell_ -= cellIsAxisPrism;
            nProtected_ -= nAxisPrims;



            // YO
            // Unprotect cells that have a refinement history
            const labelList& visibleCells =
                meshCutter_->history().visibleCells();
            forAll(visibleCells, celli)
            {
                if (visibleCells[celli] >= 0)
                {
                    if (protectedCell_.get(celli))
                    {
                        protectedCell_.unset(celli);
                        nProtected_--;
                    }
                }
            }

            // Look for former prisms cells on level 0. Since they do not have
            // a refinement history, they will be incorrectly protected.
            // In order to find former prisms, look for protectedCells with
            // exactly 2 points on the axis. If the number of refined
            // neighbours matches with the number of faces in excess of 5,
            // the cell is (most probably) a former prism.

            {
                // First, look for a wedge polyPatch, in order to get the wedge
                // plane normal and axis.
                vector axis(0,0,0);
                vector centreNormal(0,0,0);
                bool foundWedge = false;

                forAll(mesh_.boundaryMesh(), patchi)
                {
                    if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchi]))
                    {
                        foundWedge = true;

                        const wedgePolyPatch& wedgePatch =
                            refCast<const wedgePolyPatch>
                            (
                                mesh_.boundaryMesh()[patchi]
                            );

                        axis = wedgePatch.axis();
                        centreNormal = wedgePatch.centreNormal();
                        break;
                    }
                }

                if (!foundWedge)
                {
                    FatalErrorInFunction
                        << "Number of geometric dimensions and solution dimensions "
                        << "correspond to an axisymmetric case, but no wedgePolyPatch "
                        << " found in boundaryMesh." << nl
                        << exit(FatalError);
                }

                // Get radial direction on the mesh
                vector radialVector = axis ^ centreNormal;

                // Translate as a component. Since wedge meshes should be constrained
                // around one of the planes XY, XZ or YZ, we can directly get the
                // distance to the axis

                direction dir = 0;
                scalar maxComp = mag(radialVector[0]);

                if (mag(radialVector[1]) > maxComp)
                {
                    dir = 1;
                    maxComp = mag(radialVector[1]);
                }
                if (mag(radialVector[2]) > maxComp)
                {
                    dir = 2;
                }

                // Use the level0Edge length to calculate a safety margin for
                // detecting points on the axis.
                // We should be able to rely on level0, since we only check
                // unrefined cells.
                scalar tolerance = meshCutter_->level0EdgeLength() * 1e-6;

                const labelListList& cellPoints = mesh_.cellPoints();

                // Get max level across faces
                labelList maxFaceLevel(mesh_.nFaces());
                const labelList& cellLevel = meshCutter_->cellLevel();

                forAll(maxFaceLevel, facei)
                {
                    maxFaceLevel[facei] = cellLevel[mesh_.faceOwner()[facei]];
                }

                forAll(mesh_.faceNeighbour(), facei)
                {
                    maxFaceLevel[facei] =
                        max
                        (
                            maxFaceLevel[facei],
                            cellLevel[mesh_.faceNeighbour()[facei]]
                        );
                }
                syncTools::syncFaceList
                (
                    mesh_,
                    maxFaceLevel,
                    combineMaxOp<label>()
                );

                forAll(protectedCell_, celli)
                {
                    if (protectedCell_.get(celli))
                    {
                        // Check if the cell has exactly 2 points on the axis
                        label numPointsOnAxis = 0;
                        const labelList& cPoints = cellPoints[celli];
                        forAll(cPoints, pointi)
                        {
                            if
                            (
                                mag(mesh_.points()[cPoints[pointi]][dir])
                              < tolerance
                            )
                            {
                                numPointsOnAxis++;
                            }
                        }

                        if (numPointsOnAxis == 2)
                        {
                            // Calculate number of higher level neighbour cells.
                            // Use cellFaces instead of cellCells to deal with
                            // processor patches.
                            label numNeighboursWithHigherLevel = 0;

                            const labelList& cellFaces = mesh_.cells()[celli];
                            forAll(cellFaces, cellFacei)
                            {
                                if
                                (
                                    cellLevel[celli]
                                  < maxFaceLevel[cellFaces[cellFacei]]
                                )
                                {
                                    numNeighboursWithHigherLevel++;
                                }
                            }

                            // If the cell was a neighbour to refined cells, then each
                            // pair of refined neighbours introduces an additional face
                            if
                            (
                                2*(cellFaces.size() - 5)
                             == numNeighboursWithHigherLevel
                            )
                            {
                                protectedCell_.unset(celli);
                                nProtected_--;
                            }
                        }
                    }
                }
            }
        }
    }
    setProtectedCells();

    if (returnReduce(nProtected_, sumOp<label>()))
    {
        cellSet protectedCells(mesh_, "protectedCells", nProtected_);
        forAll(protectedCell_, celli)
        {
            if (protectedCell_.get(celli))
            {
                protectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected_, sumOp<label>())
            << " cells that are protected from refinement." << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshHexRefiner::~fvMeshHexRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshHexRefiner::refine
(
    const scalarField& error,
    const labelList& maxCellLevel,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar unrefineLevel
)
{
    readDict(this->dict_);
    bool hasChanged = false;

    if (preUpdate())
    {
        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(mesh_.nCells());

        if (canRefine(true))
        {
            labelList maxRefinement(maxCellLevel);
            setMaxCellLevel(maxRefinement);

            // Determine candidates for refinement (looking at field only)
            selectRefineCandidates
            (
                lowerRefineLevel,
                upperRefineLevel,
                error,
                refineCell
            );

            // Extend with a buffer layer to prevent neighbouring points
            // being unrefined.
            for (label i = 0; i < nRefinementBufferLayers_; i++)
            {
                extendMarkedCells(refineCell, maxRefinement, i == 0, true);
            }

            forAll(protectedPatches_, patchi)
            {
                const polyPatch& p =
                    mesh_.boundaryMesh()[protectedPatches_[patchi]];
                forAll(p.faceCells(), facei)
                {
                    label own = mesh_.faceOwner()[facei + p.start()];
                    refineCell.set(own, false);
                }
            }
            if (nProtected_ > 0)
            {
                forAll(protectedCell_, celli)
                {
                    if (protectedCell_.get(celli))
                    {
                        refineCell.set(celli, false);
                    }
                }
            }

            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells_,
                    maxRefinement,
                    refineCell
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                meshCutter_->history().active() = true;
                isRefining_ = true;

                // Refine/update mesh and map fields
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                // Update refineCell. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, celli)
                    {
                        label oldCelli = cellMap[celli];

                        if (oldCelli < 0)
                        {
                            newRefineCell.set(celli, 1);
                        }
                        else if (reverseCellMap[oldCelli] != celli)
                        {
                            newRefineCell.set(celli, 1);
                        }
                        else
                        {
                            newRefineCell.set(celli, refineCell.get(oldCelli));
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }

                hasChanged = true;
                isRefining_ = false;
            }
        }


        if (canUnrefine(true))
        {
            // Extend with a buffer layer to prevent neighbouring points
            // being unrefined.
            for (label i = 0; i < nUnrefinementBufferLayers_; i++)
            {
                extendMarkedCellsAcrossFaces(refineCell);
            }

            if (nProtected_ > 0)
            {
                forAll(protectedCell_, celli)
                {
                    if (protectedCell_.get(celli))
                    {
                        refineCell.set(celli, true);
                    }
                }
            }

            forAll(protectedPatches_, patchi)
            {
                const polyPatch& p =
                    mesh_.boundaryMesh()[protectedPatches_[patchi]];
                forAll(p.faceCells(), facei)
                {
                    label own = mesh_.faceOwner()[facei + p.start()];
                    refineCell.set(own, true);
                }
            }

            // Select unrefineable points that are not marked in refineCell
            labelList elemsToUnrefine
            (
                meshCutter_->selectUnrefineElems
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(error)
                )
            );

            label nSplitElems = returnReduce
            (
                elemsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitElems > 0)
            {
                isUnrefining_ = true;

                // Refine/update mesh
                unrefine(elemsToUnrefine);

                hasChanged = true;
                isUnrefining_ = false;
            }
        }

        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occasionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<hexRefRefinementHistory&>
            (
                meshCutter_->history()
            ).compact();
        }

        reduce(hasChanged, orOp<bool>());
        if (balance())
        {
            hasChanged = true;
        }

        mesh_.topoChanging(hasChanged);
        if (hasChanged)
        {
            // Reset moving flag (if any). If not using inflation we'll not
            // move, if are using inflation any follow on movePoints will set
            // it.
            mesh_.moving(false);

            // Make sure all processors have the correct instance
            mesh_.setInstance(mesh_.time().timeName());
        }
    }

    return hasChanged;
}


bool Foam::fvMeshHexRefiner::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexRef&>(meshCutter_()).setInstance(mesh_.facesInstance());

    bool writeOk =
        fvMeshRefiner::writeObject(fmt, ver, cmp, write)
     && meshCutter_->write();

    if (returnReduce(nProtected_, sumOp<label>()) > 0)
    {
        cellSet protectedCells(mesh_, "protectedCells", nProtected_);
        forAll(protectedCell_, celli)
        {
            if (protectedCell_.get(celli))
            {
                protectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected_, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << protectedCells.name()
            << "." << endl;

        protectedCells.write();
    }

    return writeOk;
}


// ************************************************************************* //
