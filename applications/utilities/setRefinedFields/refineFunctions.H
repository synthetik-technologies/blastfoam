template<class Type>
class combineMaxOp
{
    public:
    void operator()(Type& x, const Type& y) const
    {
        x = max(x, y);
    }
};


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

            error[own] = max(error[own], e);
            error[nei] = max(error[nei], e);
        }

        // Boundary faces
        forAll(error.boundaryField(), patchi)
        {
            if (error.boundaryField()[patchi].coupled())
            {
                const fvPatch& p = field.boundaryField()[patchi].patch();

                const labelUList& faceCells = p.faceCells();
                scalarField fp
                (
                    field.boundaryField()[patchi].patchInternalField()
                );

                scalarField fn
                (
                    field.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(faceCells, facei)
                {
                    scalar e =
                        pos0(mag(fp[facei] - fn[facei]) - 1e-6)
                      - 10.0*neg(mag(fp[facei] - fn[facei]) - 1e-6);
                    error[faceCells[facei]]=
                        max(error[faceCells[facei]], e);
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

void extendMaxCellLevel
(
    const fvMesh& mesh,
    labelList& cells,
    labelList& maxCellLevel,
    const label level,
    const bool top
)
{
    // Mark faces using any marked cell
    boolList markedFace(mesh.nFaces(), false);
    PackedBoolList markedCell(mesh.nCells(), false);

    forAll(cells, i)
    {
        label celli = cells[i];
        markedCell.set(celli, true);
        const cell& cFaces = mesh.cells()[celli];

        forAll(cFaces, i)
        {
            markedFace[cFaces[i]] = true;
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

    cells.resize(mesh.nCells());
    label i = 0;
    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            cells[i++] = celli;
            maxCellLevel[celli] = max(maxCellLevel[celli], level);
        }
    }
    cells.resize(i);
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


autoPtr<mapPolyMesh> unrefine
(
    fvMesh& mesh,
    const labelList& splitElems,
    PackedBoolList& protectedCell,
    hexRef& meshCutter
)
{
    polyTopoChange meshMod(mesh);

    // Play refinement commands into mesh changer.
    meshCutter.setUnrefinement(splitElems, meshMod);


    /// Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(0);
    meshCutter.calcFaceToSplitPoint(splitElems, faceToSplitPoint);


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
        labelList elemsToUnrefine
        (
            meshCutter.selectUnrefineElems
            (
                unrefineLevel,
                refineCell,
                maxCellField(mesh, vFld)
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
            unrefine
            (
                mesh,
                elemsToUnrefine,
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
    PackedBoolList& protectedCells,
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

            if (!protectedCells.get(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    nAnchors[celli]++;

                    if (nAnchors[celli] > 8)
                    {
                        protectedCells.set(celli, 1);
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
            if (nAnchors == 2)
            {
                protectedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCells.set(mesh.faceOwner()[facei], 1);
                nProtected++;
                protectedCells.set(mesh.faceNeighbour()[facei], 1);
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
                protectedCells.set(mesh.faceOwner()[facei], 1);
                nProtected++;
            }
        }

        bool isAxisym = (mesh.nGeometricD() == 2 && mesh.nSolutionD() == 3);

        if (!isAxisym)
        {
            // Also protect any cells that are less than hex
            forAll(mesh.cells(), celli)
            {
                const cell& cFaces = mesh.cells()[celli];

                if (cFaces.size() < 6)
                {
                    if (protectedCells.set(celli, 1))
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
                            if (protectedCells.set(celli, 1))
                            {
                                nProtected++;
                            }
                            break;
                        }
                    }
                }
            }

            // Check cells for 8 corner points
            checkEightAnchorPoints
            (
                mesh,
                protectedCells,
                nProtected,
                meshCutter
            );
        }
        else
        {
            PackedBoolList cellIsAxisPrism(mesh.nCells(), false);
            label nAxisPrims = 0;

            // Do not protect prisms on the axis
            forAll(mesh.cells(), celli)
            {
                const cell& cFaces = mesh.cells()[celli];

                if (cFaces.size() == 5)
                {
                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];
                        if
                        (
                            !mesh.isInternalFace(facei)
                         && mesh.faces()[facei].size() != 4
                        )
                        {
                            label patchi =
                                mesh.boundaryMesh().whichPatch(facei);
                            if
                            (
                                isA<wedgePolyPatch>(mesh.boundaryMesh()[patchi])
                            )
                            {
                                if (protectedCells.set(celli, 1))
                                {
                                    nProtected++;
                                    break;
                                }
                            }
                        }
                    }
                    if (!protectedCells.get(celli))
                    {
                        cellIsAxisPrism.set(celli, 1);
                        nAxisPrims++;
                    }
                }
                else if (cFaces.size() < 5)
                {
                    if (protectedCells.set(celli, 1))
                    {
                        nProtected++;
                    }
                }
            }

            // Check cells for 8 corner points
            checkEightAnchorPoints
            (
                mesh,
                protectedCells,
                nProtected,
                meshCutter
            );

            // Unprotect prism cells on the axis
            protectedCells -= cellIsAxisPrism;
            nProtected -= nAxisPrims;



            // YO
            // Unprotect cells that have a refinement history
            const labelList& visibleCells =
                meshCutter.history().visibleCells();
            forAll(visibleCells, celli)
            {
                if (visibleCells[celli] >= 0)
                {
                    if (protectedCells[celli])
                    {
                        protectedCells.unset(celli);
                        nProtected--;
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

                forAll(mesh.boundaryMesh(), patchi)
                {
                    if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchi]))
                    {
                        foundWedge = true;

                        const wedgePolyPatch& wedgePatch =
                            refCast<const wedgePolyPatch>
                            (
                                mesh.boundaryMesh()[patchi]
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
                vector radialVector = axis^centreNormal;

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
                scalar tolerance = meshCutter.level0EdgeLength()*1e-6;

                const labelListList& cellPoints = mesh.cellPoints();

                // Get max level across faces
                labelList maxFaceLevel(mesh.nFaces());
                const labelList& cellLevel = meshCutter.cellLevel();

                forAll(maxFaceLevel, facei)
                {
                    maxFaceLevel[facei] = cellLevel[mesh.faceOwner()[facei]];
                }

                forAll(mesh.faceNeighbour(), facei)
                {
                    maxFaceLevel[facei] = max(maxFaceLevel[facei], cellLevel[mesh.faceNeighbour()[facei]]);
                }
                syncTools::syncFaceList(mesh, maxFaceLevel, combineMaxOp<label>());

                forAll(protectedCells, celli)
                {
                    if (protectedCells[celli])
                    {
                        // Check if the cell has exactly 2 points on the axis
                        label numPointsOnAxis = 0;
                        const labelList& cPoints = cellPoints[celli];
                        forAll(cPoints, pointi)
                        {
                            if
                            (
                                mag(mesh.points()[cPoints[pointi]][dir])
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
                            label thisCellLevel = cellLevel[celli];
                            label numNeighboursWithHigherLevel = 0;

                            const labelList& cellFaces = mesh.cells()[celli];
                            forAll(cellFaces, cellFacei)
                            {
                                if
                                (
                                    thisCellLevel
                                  < maxFaceLevel[cellFaces[cellFacei]]
                                )
                                {
                                    numNeighboursWithHigherLevel++;
                                }
                            }

                            // If the cell was a neighbour to refined cells, then each
                            // pair of refined neighbours introduces an additional face
                            if (2*(cellFaces.size()-5) == numNeighboursWithHigherLevel)
                            {
                                protectedCells.unset(celli);
                                nProtected--;
                            }
                        }
                    }
                }
            }
        } //-YO

        // Protect cells that will cause a failure (from snappyHexMesh)
        boolList protectedFaces(mesh.nFaces(), false);

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
                }
            }
            if (nAnchors == 2)
            {
                protectedFaces[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh, protectedFaces, orEqOp<bool>());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (protectedFaces[facei])
            {
                protectedCells.set(mesh.faceOwner()[facei], 1);
                nProtected++;
                protectedCells.set(mesh.faceNeighbour()[facei], 1);
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
            if (protectedFaces[facei])
            {
                protectedCells.set(mesh.faceOwner()[facei], 1);
                nProtected++;
            }
        }
    }


    if (returnReduce(nProtected, sumOp<label>()) == 0)
    {
        protectedCells.clear();
    }
    else
    {

        cellSet pCells(mesh, "protectedCells", nProtected);
        forAll(protectedCells, celli)
        {
            if (protectedCells[celli])
            {
                pCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << pCells.name()
            << "." << endl;

        pCells.write();
    }
}

// ************************************************************************* //
