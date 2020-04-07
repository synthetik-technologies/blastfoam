// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool setRegion
(
    const label nCells,
    const label nOldCells,
    const label level,
    const labelList& cells,
    const labelList& cellLevel,
    volScalarField& error
)
{
    const labelUList& owner = error.mesh().owner();
    const labelUList& neighbour = error.mesh().neighbour();

    if (nCells == 0)
    {
        return false;
    }

    bool set(nCells != nOldCells);

    forAll(cells, celli)
    {
        const labelList& faces(error.mesh().cells()[cells[celli]]);
        forAll(faces, facei)
        {
            if (faces[facei] < error.mesh().nInternalFaces())
            {
                label own = owner[faces[facei]];
                label nei = neighbour[faces[facei]];

                if (cellLevel[cells[celli]] < level)
                {
                    error[own] = 0.1;
                    error[nei] = 0.1;
                }
                if (cellLevel[own] > level)
                {
                    error[own] = max(error[own], -10.0);
                }
                if (cellLevel[nei] > level)
                {
                    error[nei] = max(error[nei], -10.0);
                }
            }
        }
    }
    return set;
}


void calcFaceDiff
(
    volScalarField& error,
    const PtrList<volScalarField>& fields
)
{
    volScalarField errorOrig(error);
    const labelUList& owner = error.mesh().owner();
    const labelUList& neighbour = error.mesh().neighbour();
    const label nInternalFaces = error.mesh().nInternalFaces();

    forAll(fields, fieldi)
    {
        const volScalarField& field = fields[fieldi];
        for (label facei = 0; facei < nInternalFaces; facei++)
        {
            label own = owner[facei];
            label nei = neighbour[facei];
            scalar e =
                pos0(mag(field[own] - field[nei]) - 1e-6)
              - 10.0*neg(mag(field[own] - field[nei]) - 1e-6);

            if (errorOrig[own] > 0)
            {
                error[own] = max(error[own], e);
            }
            if (errorOrig[nei] > 0)
            {
                error[nei] = max(error[nei], e);
            }
        }

        // Boundary faces
        forAll(error.boundaryField(), patchi)
        {
            if (error.boundaryField()[patchi].coupled())
            {
                const fvPatch& p = field.boundaryField()[patchi].patch();

                const labelUList& faceCells = p.faceCells();
                const scalarField fp = field.boundaryField()[patchi].patchInternalField();

                const scalarField fn
                (
                    field.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(faceCells, facei)
                {
                    scalar e =
                        pos0(mag(fp[facei] - fn[facei]) - 1e-6)
                      - 10.0*neg(mag(fp[facei] - fn[facei]) - 1e-6);
                    if (errorOrig[faceCells[facei]] > 0)
                    {
                        error[faceCells[facei]]=
                            max(error[faceCells[facei]], e);
                    }
                }
            }
        }
    }
}


Foam::label count
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
    }

    return n;
}


void calculateProtectedCells
(
    const fvMesh& mesh,
    PackedBoolList& unrefineableCell,
    const hexRef& meshCutter,
    const PackedBoolList& protectedCell
)
{
    if (protectedCell.empty())
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter.cellLevel();

    unrefineableCell = protectedCell;

    // Get neighbouring cell level
    labelList neiLevel(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiLevel[facei-mesh.nInternalFaces()] =
            cellLevel[mesh.faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(mesh.nFaces(), false);

        forAll(mesh.faceNeighbour(), facei)
        {
            label own = mesh.faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            label nei = mesh.faceNeighbour()[facei];
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
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            label own = mesh.faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            if
            (
                ownProtected
             && (neiLevel[facei-mesh.nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh, seedFace, orEqOp<bool>());


        // Extend unrefineableCell
        bool hasExtended = false;

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = mesh.faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }

                label nei = mesh.faceNeighbour()[facei];
                if (unrefineableCell.get(nei) == 0)
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            if (seedFace[facei])
            {
                label own = mesh.faceOwner()[facei];
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
autoPtr<mapPolyMesh> refine
(
    fvMesh& mesh,
    const labelList& cellsToRefine,
    PackedBoolList& protectedCell,
    hexRef& meshCutter
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh);

    // Play refinement commands into mesh changer.
    meshCutter.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh.globalData().nTotalCells() << " cells." << endl;

    // Update fields
    mesh.updateMesh(map);

    // Update numbering of cells/vertices.
    meshCutter.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell.size())
    {
        PackedBoolList newProtectedCell(mesh.nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            newProtectedCell.set(celli, protectedCell.get(oldCelli));
        }
        protectedCell.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter.checkRefinementLevels(-1, labelList(0));

    return map;
}

scalarField maxPointField
(
    const fvMesh& mesh,
    const scalarField& pFld
)
{
    scalarField vFld(mesh.nCells(), -great);

    forAll(mesh.pointCells(), pointi)
    {
        const labelList& pCells = mesh.pointCells()[pointi];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointi]);
        }
    }
    return vFld;
}


scalarField maxCellField
(
    const fvMesh& mesh,
    const scalarField& vFld
)
{
    scalarField pFld(mesh.nPoints(), -great);

    forAll(mesh.pointCells(), pointi)
    {
        const labelList& pCells = mesh.pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = max(pFld[pointi], pFld[pCells[i]]);
        }
    }
    return pFld;
}


scalarField cellToPoint
(
    const fvMesh& mesh,
    const scalarField& vFld
)
{
    scalarField pFld(mesh.nPoints());

    forAll(mesh.pointCells(), pointi)
    {
        const labelList& pCells = mesh.pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }
    return pFld;
}


scalarField calcError
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
)
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


void selectRefineCandidates
(
    const fvMesh& mesh,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
)
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            mesh,
            calcError
            (
                cellToPoint(mesh, vFld),
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
}


labelList selectRefineCells
(
    const fvMesh& mesh,
    const label maxCells,
    const labelList& maxRefinement,
    const PackedBoolList& candidateCell,
    const hexRef& meshCutter,
    const PackedBoolList& protectedCell
)
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - mesh.globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCell;
    calculateProtectedCells
    (
        mesh,
        unrefineableCell,
        meshCutter,
        protectedCell
    );


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
        // Sort by error? For now just truncate.
        for (label level = 0; level < max(maxRefinement); level++)
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
        meshCutter.consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << mesh.globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


void extendMarkedCells
(
    const fvMesh& mesh,
    PackedBoolList& markedCell,
    const labelList& maxRefinement,
    const labelList& cellLevel,
    const bool top
)
{
    // Mark faces using any marked cell
    boolList markedFace(mesh.nFaces(), false);

    forAll(markedCell, celli)
    {
        if
        (
            markedCell.get(celli)
         && (maxRefinement[celli] > cellLevel[celli] || !top)
        )
        {
            const cell& cFaces = mesh.cells()[celli];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh.faceOwner()[facei], 1);
            markedCell.set(mesh.faceNeighbour()[facei], 1);
        }
    }
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh.faceOwner()[facei], 1);
        }
    }
}


void checkEightAnchorPoints
(
    const fvMesh& mesh,
    PackedBoolList& protectedCell,
    label& nProtected,
    const hexRef& meshCutter
)
{
    const labelList& cellLevel = meshCutter.cellLevel();
    const labelList& pointLevel = meshCutter.pointLevel();

    labelList nAnchorPoints(mesh.nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = mesh.pointCells(pointi);

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

labelList selectUnrefinePointsEdges
(
    const fvMesh& mesh,
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld,
    const hexRef& meshCutter
)
{
    // All points that can be unrefined
    const labelList splitPointsEdges(meshCutter.getSplitPointsEdges());

    DynamicList<label> newSplitPointsEdges(splitPointsEdges.size());

    forAll(splitPointsEdges, i)
    {
        if (meshCutter.useEdges())
        {
            label edgei = splitPointsEdges[i];

            const edge& e = mesh.edges()[edgei];

            forAll(e, j)
            {
                label pointi = e[j];

                bool hasMarked = true;

                // Whas is the meaning of unrefineLevel?
                // Does it makes sense to compare it with
                // the field value instead of the cell level as the name suggests?
                if (pFld[pointi] < unrefineLevel)
                {
                    hasMarked = false;

                    // Check that all cells are not marked
                    const labelList& pCells = mesh.pointCells()[pointi];



                    forAll(pCells, pCelli)
                    {
                        if (markedCell.get(pCells[pCelli]) == 1)
                        {
                            hasMarked = true;
                            break;
                        }
                    }
                }

                if (!hasMarked)
                {
                    newSplitPointsEdges.append(edgei);
                    break;
                }
            }
        }
        else
        {
            label pointi = splitPointsEdges[i];

            if (pFld[pointi] < unrefineLevel)
            {
                // Check that all cells are not marked
                const labelList& pCells = mesh.pointCells()[pointi];

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
                    newSplitPointsEdges.append(pointi);
                }
            }
        }
    }


    newSplitPointsEdges.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter.consistentUnrefinement
        (
            newSplitPointsEdges,
            false
        )
    );
    if (meshCutter.useEdges())
    {
        Info<< "Selected "
            << returnReduce(consistentSet.size(), sumOp<label>())
            << " split edges out of a possible "
            << returnReduce(splitPointsEdges.size(), sumOp<label>())
            << "." << endl;
    }
    else
    {
        Info<< "Selected "
            << returnReduce(consistentSet.size(), sumOp<label>())
            << " split points out of a possible "
            << returnReduce(splitPointsEdges.size(), sumOp<label>())
            << "." << endl;
    }

    return consistentSet;
}


autoPtr<mapPolyMesh> unrefine
(
    fvMesh& mesh,
    const labelList& splitPointsEdges,
    PackedBoolList& protectedCell,
    hexRef& meshCutter
)
{
    polyTopoChange meshMod(mesh);

    // Play refinement commands into mesh changer.
    meshCutter.setUnrefinement(splitPointsEdges, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    if (meshCutter.useEdges())
    {
        Map<label> faceToSplitPoint(2*splitPointsEdges.size());

        forAll(splitPointsEdges, i)
        {
            label edgei = splitPointsEdges[i];

            const edge& e = mesh.edges()[edgei];

            forAll(e, j)
            {
                label pointi = e[j];

                const labelList& pFaces = mesh.pointFaces()[pointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], pointi);
                }
            }
        }
    }
    else
    {
        Map<label> faceToSplitPoint(3*splitPointsEdges.size());

        forAll(splitPointsEdges, i)
        {
            label pointi = splitPointsEdges[i];

            const labelList& pEdges = mesh.pointEdges()[pointi];

            forAll(pEdges, j)
            {
                label otherPointi =
                    mesh.edges()[pEdges[j]].otherVertex(pointi);

                const labelList& pFaces = mesh.pointFaces()[otherPointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], otherPointi);
                }
            }
        }
    }


    // Change mesh and generate map.
    // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh.globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    mesh.updateMesh(map);


    // Update numbering of cells/vertices.
    meshCutter.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell.size())
    {
        PackedBoolList newProtectedCell(mesh.nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            if (oldCelli >= 0)
            {
                newProtectedCell.set(celli, protectedCell.get(oldCelli));
            }
        }
        protectedCell.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter.checkRefinementLevels(-1, labelList(0));

    return map;
}


bool update
(
    fvMesh& mesh,
    const volScalarField& vFld,
    const dictionary& refineDict,
    hexRef& meshCutter,
    PackedBoolList& protectedCell,
    const labelList& maxRefinement
)
{
    bool hasChanged = false;

    label maxCells(refineDict.lookupOrDefault<label>("maxCells", labelMax));

    if (maxCells <= 0)
    {
        FatalErrorInFunction
            << "Illegal maximum number of cells " << maxCells << nl
            << "The maxCells setting in the setFieldsDict should"
            << " be > 0." << nl
            << exit(FatalError);
    }


    if (min(maxRefinement) <= 0)
    {
        FatalErrorInFunction
            << "Illegal maximum refinement level " << min(maxRefinement) << nl
            << "The maxCells setting in the setFields should"
            << " be > 0." << nl
            << exit(FatalError);
    }

    const scalar lowerRefineLevel = 0.5;
    const scalar unrefineLevel = -0.5;

    const label nBufferLayers =
        readLabel(refineDict.lookup("nBufferLayers"));

    // Cells marked for refinement or otherwise protected from unrefinement.
    PackedBoolList refineCell(mesh.nCells());

    // Determine candidates for refinement (looking at field only)
    selectRefineCandidates
    (
        mesh,
        lowerRefineLevel,
        great,
        vFld,
        refineCell
    );

    if (mesh.globalData().nTotalCells() < maxCells)
    {
        // Extend with a buffer layer to prevent neighbouring points
        // being unrefined.
        for (label i = 0; i < nBufferLayers; i++)
        {
            extendMarkedCells
            (
                mesh,
                refineCell,
                maxRefinement,
                meshCutter.cellLevel(),
                i == 0
            );
        }

        // Select subset of candidates. Take into account max allowable
        // cells, refinement level, protected cells.
        labelList cellsToRefine
        (
            selectRefineCells
            (
                mesh,
                maxCells,
                maxRefinement,
                refineCell,
                meshCutter,
                protectedCell
            )
        );

        label nCellsToRefine = returnReduce
        (
            cellsToRefine.size(), sumOp<label>()
        );

        if (nCellsToRefine > 0)
        {
            // Refine/update mesh and map fields
            autoPtr<mapPolyMesh> map =
                refine
                (
                    mesh,
                    cellsToRefine,
                    protectedCell,
                    meshCutter
                );

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

    {
        // Select unrefineable points that are not marked in refineCell
        labelList pointsEdgesToUnrefine
        (
            selectUnrefinePointsEdges
            (
                mesh,
                unrefineLevel,
                refineCell,
                maxCellField(mesh, vFld),
                meshCutter
            )
        );

        label nSplitPointsEdges = returnReduce
        (
            pointsEdgesToUnrefine.size(),
            sumOp<label>()
        );

        if (nSplitPointsEdges > 0)
        {
            // Refine/update mesh
            unrefine
            (
                mesh,
                pointsEdgesToUnrefine,
                protectedCell,
                meshCutter
            );

            hasChanged = true;
        }
    }

    return hasChanged;
}

void initialize
(
    fvMesh& mesh,
    PackedBoolList& protectedCell,
    hexRef& meshCutter
)
{
    const labelList& cellLevel = meshCutter.cellLevel();
    const labelList& pointLevel = meshCutter.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(mesh.nCells(), 0);

    label nProtected = 0;

    forAll(mesh.pointCells(), pointi)
    {
        const labelList& pCells = mesh.pointCells()[pointi];

        forAll(pCells, i)
        {
            label celli = pCells[i];

            if (!protectedCell.get(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    nAnchors[celli]++;

                    if (nAnchors[celli] > 8)
                    {
                        protectedCell.set(celli, 1);
                        nProtected++;
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
        labelList neiLevel(mesh.nFaces());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[mesh.faceNeighbour()[facei]];
        }
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            neiLevel[facei] = cellLevel[mesh.faceOwner()[facei]];
        }
        syncTools::swapFaceList(mesh, neiLevel);


        boolList protectedFace(mesh.nFaces(), false);

        forAll(mesh.faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[mesh.faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = mesh.faces()[facei];

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

        syncTools::syncFaceList(mesh, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell.set(mesh.faceOwner()[facei], 1);
                nProtected++;
                protectedCell.set(mesh.faceNeighbour()[facei], 1);
                nProtected++;
            }
        }
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            if (protectedFace[facei])
            {
                protectedCell.set(mesh.faceOwner()[facei], 1);
                nProtected++;
            }
        }

        // Also protect any cells that are less than hex
        forAll(mesh.cells(), celli)
        {
            const cell& cFaces = mesh.cells()[celli];

            if (cFaces.size() < 6)
            {
                if (protectedCell.set(celli, 1))
                {
                    nProtected++;
                }
            }
            else
            {
                forAll(cFaces, cFacei)
                {
                    if (mesh.faces()[cFaces[cFacei]].size() < 4)
                    {
                        if (protectedCell.set(celli, 1))
                        {
                            nProtected++;
                        }
                        break;
                    }
                }
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(mesh, protectedCell, nProtected, meshCutter);
    }

    if (returnReduce(nProtected, sumOp<label>()) == 0)
    {
        protectedCell.clear();
    }
    else
    {

        cellSet ProtectedCells(mesh, "protectedCells", nProtected);
        forAll(protectedCell, celli)
        {
            if (protectedCell[celli])
            {
                ProtectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << ProtectedCells.name()
            << "." << endl;

        ProtectedCells.write();
    }
}

// ************************************************************************* //
