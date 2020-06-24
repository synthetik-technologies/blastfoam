/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020  Jeff Heylmun:   Modified original dynamicRefineBalanceFvMesh class
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

#include "adaptiveFvMesh.H"
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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adaptiveFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, adaptiveFvMesh, IOobject);

    // Helper class for accessing max cell level of faces accross processor patches
    template<class Type>
    class combineMaxOp
    {
        public:
        void operator()(Type& x, const Type& y) const
        {
            x = max(x, y);
        }
    };
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::adaptiveFvMesh::topParentID(const label p) const
{
    if (p >= meshCutter().history().splitCells().size())
    {
        return p;
    }
    label nextP = meshCutter().history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::adaptiveFvMesh::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // debug also serves to get-around Clang compiler trying to optimsie
        // out this forAll loop under O3 optimisation
        if (debug)
        {
            Info<< "n=" << n << endl;
        }
    }

    return n;
}


void Foam::adaptiveFvMesh::calculateProtectedCells
(
    PackedBoolList& unrefineableCell
) const
{
    if (protectedCell_.empty())
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_->cellLevel();

    unrefineableCell = protectedCell_;

    // Get neighbouring cell level
    labelList neiLevel(nFaces()-nInternalFaces());

    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    {
        neiLevel[facei-nInternalFaces()] = cellLevel[faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(*this, neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(nFaces(), false);

        forAll(faceNeighbour(), facei)
        {
            label own = faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            label nei = faceNeighbour()[facei];
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
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            label own = faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            if
            (
                ownProtected
             && (neiLevel[facei-nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(*this, seedFace, orEqOp<bool>());


        // Extend unrefineableCell
        bool hasExtended = false;

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }

                label nei = faceNeighbour()[facei];
                if (unrefineableCell.get(nei) == 0)
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = faceOwner()[facei];
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


void Foam::adaptiveFvMesh::readDict()
{
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );
    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));

    if (!refineDict.found("correctFluxes"))
    {
        return;
    }

    List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        refineDict.lookup("correctFluxes")
    );
    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::adaptiveFvMesh::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_->setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            label oldFacei = map().faceMap()[facei];

            if (oldFacei >= nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << faceCentres()[facei]
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
    this->updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, facei)
        {
            label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                label masterFacei = reverseFaceMap[oldFacei];

                if (masterFacei < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << facei << abort(FatalError);
                }
                else if (masterFacei != facei)
                {
                    masterFaces.insert(masterFacei);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << iter.key()
                    << " to NaN" << endl;

                surfaceScalarField& phi = *iter();

                sigFpe::fillNan(phi.primitiveFieldRef());

                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    phi[facei] = phiU[facei];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    phi[facei] = phiU[facei];
                }
            }

            // Recalculate new boundary faces.
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();
            forAll(phiBf, patchi)
            {
                fvsPatchScalarField& patchPhi = phiBf[patchi];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchi];

                label facei = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFacei = faceMap[facei];

                    if (oldFacei == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    facei++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label facei = iter.key();

                if (isInternalFace(facei))
                {
                    phi[facei] = phiU[facei];
                }
                else
                {
                    label patchi = boundaryMesh().whichPatch(facei);

                    if (!isA<emptyPolyPatch>(boundaryMesh()[patchi]))
                    {
                        label i = facei - boundaryMesh()[patchi].start();

                        const fvsPatchScalarField& patchPhiU =
                            phiU.boundaryField()[patchi];

                        fvsPatchScalarField& patchPhi = phiBf[patchi];

                        patchPhi[i] = patchPhiU[i];
                    }
                }
            }
        }
    }

    //- map non-flux surface<Type>Fields for newly created internal faces
    const labelList& faceMap = map().faceMap();
    mapNewInternalFaces<scalar>( faceMap );
    mapNewInternalFaces<vector>( faceMap );

    // Update numbering of cells/vertices.
    meshCutter_->updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

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
Foam::adaptiveFvMesh::unrefine
(
    const labelList& splitElems
)
{
    polyTopoChange meshMod(*this);

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
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    this->updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );


            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFacei = iter.key();
                label oldPointi = iter();

                if (reversePointMap[oldPointi] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    label facei = reverseFaceMap[oldFacei];

                    if (facei >= 0)
                    {
                        if (isInternalFace(facei))
                        {
                            phi[facei] = phiU[facei];
                        }
                        else
                        {
                            label patchi = boundaryMesh().whichPatch(facei);

                            if (!isA<emptyPolyPatch>(boundaryMesh()[patchi]))
                            {
                                label i = facei - boundaryMesh()[patchi].start();

                                const fvsPatchScalarField& patchPhiU =
                                    phiU.boundaryField()[patchi];
                                fvsPatchScalarField& patchPhi = phiBf[patchi];
                                patchPhi[i] = patchPhiU[i];
                            }
                        }
                    }
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_->updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

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


Foam::scalarField
Foam::adaptiveFvMesh::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(nCells(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointi]);
        }
    }
    return vFld;
}

Foam::scalarField
Foam::adaptiveFvMesh::maxCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = max(pFld[pointi], vFld[pCells[i]]);
        }
    }
    return pFld;
}

/*
Foam::scalarField
Foam::adaptiveFvMesh::minCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = max(pFld[pointi], vFld[pCells[i]]);
        }
    }
    return pFld;
}
*/

// Simple (non-parallel) interpolation by averaging.
Foam::scalarField
Foam::adaptiveFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }
    return pFld;
}


Foam::scalarField Foam::adaptiveFvMesh::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }
    return c;
}


void Foam::adaptiveFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
////////////////////////////////////////////////////////////////////////
// Batzdorf
/*
    //-
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if (cellError[celli] > 0)
        {
            candidateCell.set(celli, 1);
        }
    }
*/
    //+
    forAll(vFld, cellI)
    {
        if ((vFld[cellI] >= lowerRefineLevel) && (vFld[cellI] <= upperRefineLevel))
        {
            candidateCell.set(cellI, 1);
        }
    }
}


Foam::labelList Foam::adaptiveFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCell
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

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
                cellLevel[celli] < maxRefinement
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
        // Sort by error? For now just truncate.
        for (label level = 0; level < maxRefinement; level++)
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
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}

// YO- This is here only to preserve compatibility with the official release.
//     It is not used by the refinement procedure, but some utilities such as
//     decomposePar rely on it.
Foam::labelList Foam::adaptiveFvMesh::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_->getSplitElems());

    DynamicList<label> newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        label pointi = splitPoints[i];

        if (pFld[pointi] < unrefineLevel)
        {
            // Check that all cells are not marked
            const labelList& pCells = pointCells()[pointi];

            bool hasMarked = false;

            forAll(pCells, pCelli)
            {
                if (markedCell.get(pCells[pCelli]))
                {
                    hasMarked = true;
                    break;
                }
            }

            if (!hasMarked)
            {
                newSplitPoints.append(pointi);
            }
        }
    }


    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_->consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}
//-YO

void Foam::adaptiveFvMesh::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(nFaces(), false);

    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            const cell& cFaces = cells()[celli];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(faceOwner()[facei], 1);
            markedCell.set(faceNeighbour()[facei], 1);
        }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(faceOwner()[facei], 1);
        }
    }
}


void Foam::adaptiveFvMesh::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    labelList nAnchorPoints(nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = pointCells(pointi);

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

                if (!protectedCell[celli])
                {
                    nAnchorPoints[celli]++;
                }
            }
        }
    }


    forAll(protectedCell, celli)
    {
        if (!protectedCell[celli] && nAnchorPoints[celli] != 8)
        {
            protectedCell.set(celli, true);
            nProtected++;
        }
    }
}


void Foam::adaptiveFvMesh::setProtectedCells()
{
    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(nCells(), 0);

    nProtected_ = 0;

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

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

    const dictionary& refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );
    wordList protectedPatches
    (
        refineDict.lookupOrDefault("protectedPatches", wordList())
    );
    if (protectedPatches.size())
    {
        forAll(protectedPatches, patchi)
        {
            const polyPatch& p = this->boundaryMesh()[protectedPatches[patchi]];
            forAll(p.faceCells(), facei)
            {
                protectedCell_.set(p.faceCells()[facei], 1);
                nProtected_++;
            }
        }
    }

    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(nFaces());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[faceNeighbour()[facei]];
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[faceOwner()[facei]];
        }
        syncTools::swapFaceList(*this, neiLevel);


        boolList protectedFace(nFaces(), false);

        forAll(faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = faces()[facei];

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

        syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(faceOwner()[facei], 1);
                nProtected_++;
                protectedCell_.set(faceNeighbour()[facei], 1);
                nProtected_++;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(faceOwner()[facei], 1);
                nProtected_++;
            }
        }

        bool isAxisym = (nGeometricD() == 2 && nSolutionD() == 3);

        if (!isAxisym)
        {
            // Also protect any cells that are less than hex
            forAll(cells(), celli)
            {
                const cell& cFaces = cells()[celli];

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
                        if (faces()[cFaces[cFacei]].size() < 4)
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
            PackedBoolList cellIsAxisPrism(nCells(), false);
            label nAxisPrims = 0;

            // Do not protect prisms on the axis
            forAll(cells(), celli)
            {
                const cell& cFaces = cells()[celli];

                if (cFaces.size() == 5)
                {
                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];
                        if
                        (
                            !isInternalFace(facei)
                        && faces()[facei].size() != 4
                        )
                        {
                            label patchi = boundaryMesh().whichPatch(facei);
                            if (isA<wedgePolyPatch>(boundaryMesh()[patchi]))
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
            const labelList& visibleCells = meshCutter_->history().visibleCells();
            forAll(visibleCells, celli)
            {
                if (visibleCells[celli] >= 0)
                {
                    if (protectedCell_[celli])
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

                forAll(boundaryMesh(), patchi)
                {
                    if (isA<wedgePolyPatch>(boundaryMesh()[patchi]))
                    {
                        foundWedge = true;

                        const wedgePolyPatch& wedgePatch =
                            refCast<const wedgePolyPatch>
                            (
                                boundaryMesh()[patchi]
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

                const labelListList& cellPoints = primitiveMesh::cellPoints();

                // Get max level across faces
                labelList maxFaceLevel(nFaces());
                const labelList& cellLevel = meshCutter_->cellLevel();

                forAll(maxFaceLevel, facei)
                {
                    maxFaceLevel[facei] = cellLevel[faceOwner()[facei]];
                }

                forAll(faceNeighbour(), facei)
                {
                    maxFaceLevel[facei] = max(maxFaceLevel[facei], cellLevel[faceNeighbour()[facei]]);
                }
                syncTools::syncFaceList(*this, maxFaceLevel, combineMaxOp<label>());

                forAll(protectedCell_, celli)
                {
                    if (protectedCell_[celli])
                    {
                        // Check if the cell has exactly 2 points on the axis
                        label numPointsOnAxis = 0;
                        const labelList& cPoints = cellPoints[celli];
                        forAll(cPoints, pointi)
                        {
                            if (mag(points()[cPoints[pointi]][dir]) < tolerance)
                            {
                                numPointsOnAxis++;
                            }
                        }

                        if (numPointsOnAxis == 2)
                        {
                            // Calculate number of higher level neighbour cells.
                            // Use cellFaces instead of cellCells to deal with
                            // processor patches.
                            label thisCellLevel = cellLevel[celli];
                            label numNeighboursWithHigherLevel = 0;

                            const labelList& cellFaces = cells()[celli];
                            forAll(cellFaces, cellFacei)
                            {
                                if (thisCellLevel < maxFaceLevel[cellFaces[cellFacei]])
                                {
                                    numNeighboursWithHigherLevel++;
                                }
                            }

                            // If the cell was a neighbour to refined cells, then each
                            // pair of refined neighbours introduces an additional face
                            if (2*(cellFaces.size()-5) == numNeighboursWithHigherLevel)
                            {
                                protectedCell_.unset(celli);
                                nProtected_--;
                            }
                        }
                    }
                }
            }
        }

        // Protect cells that will cause a failure (from snappyHexMesh)
        boolList protectedFaces(nFaces(), false);

        forAll(faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = faces()[facei];

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

        syncTools::syncFaceList(*this, protectedFaces, orEqOp<bool>());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (protectedFaces[facei])
            {
                protectedCell_.set(faceOwner()[facei], 1);
                nProtected_++;
                protectedCell_.set(faceNeighbour()[facei], 1);
                nProtected_++;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (protectedFaces[facei])
            {
                protectedCell_.set(faceOwner()[facei], 1);
                nProtected_++;
            }
        }

        //-YO
    }

    if (returnReduce(nProtected_, sumOp<label>()) == 0)
    {
        protectedCell_.clear();
    }
}


template<class T>
void Foam::adaptiveFvMesh::mapNewInternalFaces
(
    const labelList& faceMap
)
{
    return;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveFvMesh::adaptiveFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    error_(errorEstimator::New(*this, dynamicMeshDict())),
    meshCutter_(hexRef::New(*this)),
    dumpLevel_(false),
    nRefinementIterations_(0),
    nProtected_(0),
    protectedCell_(nCells(), 0),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                time().system(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    )
{
    // Read static part of dictionary
    readDict();
    setProtectedCells();

    if (returnReduce(nProtected_, sumOp<label>()))
    {
        cellSet protectedCells(*this, "protectedCells", nProtected_);
        forAll(protectedCell_, celli)
        {
            if (protectedCell_[celli])
            {
                protectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected_, sumOp<label>())
            << " cells that are protected from refinement." << endl;
    }

    //- 2D refinment does not currently work
    //  Refinement history is not compatable with the current method of
    //  tracking parent cells
    const dictionary& balanceDict =
        dynamicMeshDict().optionalSubDict("loadBalance");
    balance_ = balanceDict.lookupOrDefault("balance", true);

    if (!Pstream::parRun())
    {
        balance_ = false;
    }

    if (balance_)
    {
        // Change decomposition method if entry is present
        if (balanceDict.found("method"))
        {
            decompositionDict_.set
            (
                "method",
                balanceDict.lookup("method")
            );
        }

        // Add refinementHistory constraint
        {
            dictionary refinementHistoryDict;
            refinementHistoryDict.add("type", "refinementHistory");
            dictionary constraintsDict;
            constraintsDict.add("refinementHistory", refinementHistoryDict);
            decompositionDict_.add("constraints", constraintsDict);
        }

        decomposer_ = decompositionMethod::New(decompositionDict_);
        if (!decomposer_->parallelAware())
        {
            FatalErrorInFunction
                << "You have selected decomposition method "
                << decomposer_->typeName
                << " which is not parallel aware." << endl
                << "Please select one that is (hierarchical, ptscotch)"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adaptiveFvMesh::~adaptiveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adaptiveFvMesh::mapFields(const mapPolyMesh& mpm)
{
// DebugVar(mpm.nOldCells());
    dynamicFvMesh::mapFields(mpm);

    // Correct surface fields on introduced internal faces. These get
    // created out-of-nothing so get an interpolated value.
    mapNewInternalFaces<scalar>(mpm.faceMap());
    mapNewInternalFaces<vector>(mpm.faceMap());
    mapNewInternalFaces<sphericalTensor>(mpm.faceMap());
    mapNewInternalFaces<symmTensor>(mpm.faceMap());
    mapNewInternalFaces<tensor>(mpm.faceMap());
}


bool Foam::adaptiveFvMesh::update()
{
    //- Update error field
    error_->update();

    // Re-read dictionary. Chosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    const dictionary& refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    label refineInterval = readLabel(refineDict.lookup("refineInterval"));
    scalar beginUnrefine = refineDict.lookupOrDefault("beginUnrefine", 0.0);

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorInFunction
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }

    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if
    (
        time().timeIndex() > 0
     && time().timeIndex() % refineInterval == 0
    )
    {
        if (returnReduce(nProtected_, sumOp<label>()) > 0 && balance_)
        {
            setProtectedCells();
        }

        label maxCells
        (
            refineDict.lookupOrDefault("maxCells", labelMax)
        );

        if (maxCells <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        label maxRefinement = readLabel(refineDict.lookup("maxRefinement"));

        if (maxRefinement == 0)
        {
            topoChanging(hasChanged);

            return false;
        }
        else if (maxRefinement < 0)
        {
            FatalErrorInFunction
                << "Illegal maximum refinement level " << maxRefinement << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        const scalar lowerRefineLevel =
            readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel =
            refineDict.lookupOrDefault("upperRefineLevel", great);
        const scalar unrefineLevel = refineDict.lookupType<scalar>
        (
            "unrefineLevel"
        );
        const label nBufferLayers =
            readLabel(refineDict.lookup("nBufferLayers"));

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(nCells());

        // Determine candidates for refinement (looking at field only)
        selectRefineCandidates
        (
            lowerRefineLevel,
            upperRefineLevel,
            error_(),
            refineCell
        );

        if (globalData().nTotalCells() < maxCells)
        {
            // Extend with a buffer layer to prevent neighbouring points
            // being unrefined.
            for (label i = 0; i < nBufferLayers; i++)
            {
                extendMarkedCells(refineCell);
            }

            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells,
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
            }
        }


        if (time().value() > beginUnrefine)
        {
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
            // Select unrefineable points that are not marked in refineCell
            labelList elemsToUnrefine
            (
                meshCutter_->selectUnrefineElems
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(error_())
                )
            );

            label nSplitElems = returnReduce
            (
                elemsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitElems > 0)
            {
                // Refine/update mesh
                unrefine(elemsToUnrefine);

                hasChanged = true;
            }
        }

        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occasionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<refinementHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;

        topoChanging(hasChanged);
        if (hasChanged)
        {
            // Reset moving flag (if any). If not using inflation we'll not move,
            // if are using inflation any follow on movePoints will set it.
            moving(false);
        }
    }
    if (hasChanged)
    {
        balance();
    }
    return hasChanged;
}


void Foam::adaptiveFvMesh::balance()
{
    //Part 1 - Call normal update from dynamicRefineFvMesh
    const dictionary& balanceDict
    (
        dynamicMeshDict().optionalSubDict("loadBalance")
    );

    label balanceInterval =
        balanceDict.lookupOrDefault("balanceInterval", 1);

    if(!balance_)
    {
        return;
    }
    if
    (
        (nRefinementIterations_ % balanceInterval) != 0
     && (nRefinementIterations_ != 1)
    )
    {
        return;
    }

    //Correct values on all coupled patches
    correctBoundaries<scalar>();
    correctBoundaries<vector>();
    correctBoundaries<sphericalTensor>();
    correctBoundaries<symmTensor>();
    correctBoundaries<tensor>();

    // Part 2 - Load Balancing
    {
        const scalar allowableImbalance =
            balanceDict.lookupOrDefault
            (
                "allowableImbalance",
                0.2
            );

        //First determine current level of imbalance - do this for all
        // parallel runs with a changing mesh, even if balancing is disabled
        label nGlobalCells = globalData().nTotalCells();
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/idealNCells;

        Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;

        //If imbalanced, construct weighted coarse graph (level 0) with node
        // weights equal to their number of subcells. This partitioning works
        // as long as the number of level 0 cells is several times greater than
        // the number of processors.
        if( maxImbalance > allowableImbalance)
        {
            Info << "\n**Solver hold for redistribution at time = "  << time().timeName() << " s" << endl;

            const labelIOList& cellLevel = meshCutter().cellLevel();
            Map<label> coarseIDmap(100);
            labelList uniqueIndex(nCells(),0);

            label nCoarse = 0;

            forAll(cells(), cellI)
            {
                if
                (
                    cellLevel[cellI] > 0
                 && meshCutter().history().visibleCells()[cellI] >= 0
                )
                {
                    //YO- 2D refinement uses fixed lists with unset parents;
                    //    we need to check that the parentIndex is set
                    label parentI = meshCutter().history().parentIndex(cellI);

                    if (parentI >= 0)
                    {
                        uniqueIndex[cellI] = nCells() + topParentID
                        (
                            meshCutter().history().parentIndex(cellI)
                        );
                    }
                    else
                    {
                        uniqueIndex[cellI] = cellI;
                    }
                    //-YO
                }
                else
                {
                    uniqueIndex[cellI] = cellI;
                }

                if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
                {
                    ++nCoarse;
                }
            }

            // Convert to local sequential indexing and calculate coarse
            // points and weights
            labelList localIndex(nCells(),0);
            pointField coarsePoints(nCoarse,vector::zero);
            scalarField coarseWeights(nCoarse,0.0);
            label nRefinementDimensions(nGeometricD());

            forAll(uniqueIndex, cellI)
            {
                localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];

                // If 2D refinement (quadtree) is ever implemented, this '3'
                // should be set in general as the number of refinement
                // dimensions.
                label w = (1 << (nRefinementDimensions*cellLevel[cellI]));

                coarseWeights[localIndex[cellI]] += 1.0;
                coarsePoints[localIndex[cellI]] += C()[cellI]/w;
            }

            labelList finalDecomp = decomposer_().decompose
            (
                *this,
                localIndex,
                coarsePoints,
                coarseWeights
            );

            scalar tolDim = globalMeshData::matchTol_*bounds().mag();

            Info<< "Distributing the mesh ..." << endl;
            fvMeshDistribute distributor(*this, tolDim);

            Info<< "Mapping the fields ..." << endl;
            autoPtr<mapDistributePolyMesh> map =
                distributor.distribute(finalDecomp);

            Info<< "Distribute the map ..." << endl;
            meshCutter_->distribute(map);


            Info << "Successfully distributed mesh" << endl;

            scalarList procLoadNew (Pstream::nProcs(), 0.0);
            procLoadNew[Pstream::myProcNo()] = this->nCells();

            reduce(procLoadNew, sumOp<List<scalar> >());

            scalar overallLoadNew = sum(procLoadNew);
            scalar averageLoadNew = overallLoadNew/double(Pstream::nProcs());

            Info << "Max deviation: " << max(Foam::mag(procLoadNew-averageLoadNew)/averageLoadNew)*100.0 << " %" << endl;
        }
    }

    //Correct values on all coupled patches
    correctBoundaries<scalar>();
    correctBoundaries<vector>();
    correctBoundaries<sphericalTensor>();
    correctBoundaries<symmTensor>();
    correctBoundaries<tensor>();

    return;
}


bool Foam::adaptiveFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexRef&>(meshCutter_()).setInstance(time().timeName());

    bool writeOk =
    (
        dynamicFvMesh::writeObject(fmt, ver, cmp, write)
     && meshCutter_->write()
    );

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar(dimless, 0)
        );

        const labelList& cellLevel = meshCutter_->cellLevel();

        forAll(cellLevel, celli)
        {
            scalarCellLevel[celli] = cellLevel[celli];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }
    if (returnReduce(nProtected_, sumOp<label>()) > 0 && balance_)
    {
        cellSet protectedCells(*this, "protectedCells", nProtected_);
        forAll(protectedCell_, celli)
        {
            if (protectedCell_[celli])
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
