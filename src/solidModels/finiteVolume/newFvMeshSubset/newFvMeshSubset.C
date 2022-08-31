/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

    MJ 23/03/05 on coupled faces change the patch of the face to the
    oldInternalFaces patch.

\*---------------------------------------------------------------------------*/

#include "newFvMeshSubset.H"
#include "boolList.H"
#include "Pstream.H"
#include "emptyPolyPatch.H"
#include "demandDrivenData.H"
#include "cyclicPolyPatch.H"
#if FOAMEXTEND > 40
    #include "oversetPolyPatch.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::newFvMeshSubset::checkCellSubset() const
{
    if (!newFvMeshSubsetPtr_)
    {
        FatalErrorIn("bool newFvMeshSubset::checkCellSubset() const")
            << "Mesh subset not set.  Please set the cell map using "
            << "void setCellSubset(const labelHashSet& cellsToSubset)" << endl
            << "before attempting to access subset data"
            << abort(FatalError);

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::newFvMeshSubset::markPoints
(
    const labelList& curPoints,
    Map<label>& pointMap
)
{
    forAll (curPoints, pointI)
    {
        // Note: insert will only insert if not yet there.
        pointMap.insert(curPoints[pointI], 0);
    }
}


void Foam::newFvMeshSubset::markPoints
(
    const labelList& curPoints,
    labelList& pointMap
)
{
    forAll (curPoints, pointI)
    {
        pointMap[curPoints[pointI]] = 0;
    }
}


// Synchronize nCellsUsingFace on both sides of coupled patches. Marks
// faces that become 'uncoupled' with 3.
void Foam::newFvMeshSubset::doCoupledPatches
(
    const bool syncPar,
    labelList& nCellsUsingFace
) const
{
    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();

    label nUncoupled = 0;

    if (syncPar && Pstream::parRun())
    {
        // Send face usage across processor patches
        forAll (oldPatches, oldPatchI)
        {
            const polyPatch& pp = oldPatches[oldPatchI];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                OPstream toNeighbour
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );

                toNeighbour
                    << SubList<label>(nCellsUsingFace, pp.size(), pp.start());
            }
        }


        // Receive face usage count and check for faces that become uncoupled.
        forAll (oldPatches, oldPatchI)
        {
            const polyPatch& pp = oldPatches[oldPatchI];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                IPstream fromNeighbour
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );

                labelList nbrCellsUsingFace(fromNeighbour);

                // Combine with this side.

                forAll (pp, i)
                {
                    if
                    (
                        nCellsUsingFace[pp.start()+i] == 1
                     && nbrCellsUsingFace[i] == 0
                    )
                    {
                        // Face's neighbour is no longer there. Mark face off
                        // as coupled
                        nCellsUsingFace[pp.start()+i] = 3;
                        nUncoupled++;
                    }
                }
            }
        }
    }

    // Do same for cyclics.
    forAll (oldPatches, oldPatchI)
    {
        const polyPatch& pp = oldPatches[oldPatchI];

        if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            forAll (cycPatch, i)
            {
                label thisFaceI = cycPatch.start() + i;
                label otherFaceI = cycPatch.transformGlobalFace(thisFaceI);

                if
                (
                    nCellsUsingFace[thisFaceI] == 1
                 && nCellsUsingFace[otherFaceI] == 0
                )
                {
                    nCellsUsingFace[thisFaceI] = 3;
                    nUncoupled++;
                }
            }
        }
    }

    if (syncPar)
    {
        reduce(nUncoupled, sumOp<label>());
    }

    if (nUncoupled > 0)
    {
        Info<< "Uncoupled " << nUncoupled << " faces on coupled patches. "
            << "(processorPolyPatch, cyclicPolyPatch)" << nl
            << "You might need to run couplePatches to restore the patch face"
            << " ordering." << nl << endl;
    }
}


labelList Foam::newFvMeshSubset::subset
(
    const label nElems,
    const labelList& selectedElements,
    const labelList& subsetMap
)
{
    // Mark selected elements.
    boolList selected(nElems, false);
    forAll (selectedElements, i)
    {
        selected[selectedElements[i]] = true;
    }

    // Count subset of selected elements
    label n = 0;
    forAll (subsetMap, i)
    {
        if (selected[subsetMap[i]])
        {
            n++;
        }
    }

    // Collect selected elements
    labelList subsettedElements(n);
    n = 0;

    forAll(subsetMap, i)
    {
        if (selected[subsetMap[i]])
        {
            subsettedElements[n++] = i;
        }
    }

    return subsettedElements;
}


void Foam::newFvMeshSubset::subsetZones()
{
    // Keep all zones, even if zero size.

    const pointZoneMesh& pointZones = baseMesh().pointZones();

    // PointZones
    List<pointZone*> pZonePtrs(pointZones.size());

    forAll(pointZones, i)
    {
        const pointZone& pz = pointZones[i];

        pZonePtrs[i] = new pointZone
        (
            pz.name(),
            subset(baseMesh().allPoints().size(), pz, pointMap()),
            i,
            newFvMeshSubsetPtr_->pointZones()
        );
    }


    // FaceZones

    const faceZoneMesh& faceZones = baseMesh().faceZones();


    // Do we need to remove zones where the side we're interested in
    // no longer exists? Guess not.
    List<faceZone*> fZonePtrs(faceZones.size());

    forAll(faceZones, i)
    {
        const faceZone& fz = faceZones[i];

        // Create list of mesh faces part of the new zone
        labelList subAddressing
        (
            subset
            (
                baseMesh().allFaces().size(),
                fz,
                faceMap()
            )
        );

        // Flipmap for all mesh faces
        boolList fullFlipStatus(baseMesh().allFaces().size(), false);
        forAll(fz, j)
        {
            fullFlipStatus[fz[j]] = fz.flipMap()[j];
        }
        // Extract sub part
        boolList subFlipStatus(subAddressing.size(), false);
        forAll(subAddressing, j)
        {
            subFlipStatus[j] = fullFlipStatus[faceMap()[subAddressing[j]]];
        }

        fZonePtrs[i] = new faceZone
        (
            fz.name(),
            subAddressing,
            subFlipStatus,
            i,
            newFvMeshSubsetPtr_->faceZones()
        );
    }


    const cellZoneMesh& cellZones = baseMesh().cellZones();

    List<cellZone*> cZonePtrs(cellZones.size());

    forAll(cellZones, i)
    {
        const cellZone& cz = cellZones[i];

        cZonePtrs[i] = new cellZone
        (
            cz.name(),
            subset(baseMesh().nCells(), cz, cellMap()),
            i,
            newFvMeshSubsetPtr_->cellZones()
        );
    }


    // Add the zones
    newFvMeshSubsetPtr_->addZones(pZonePtrs, fZonePtrs, cZonePtrs);
}


void Foam::newFvMeshSubset::makeFvDictionaries()
{
    // Try accessing dictionaries
    objectRegistry obr
    (
        IOobject
        (
            this->name() + "Subset",
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary fvSchemes
    (
        IOobject
        (
            "fvSchemes",
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    if (!fvSchemes.headerOk())
    {
        Info << "Cannot read " << fvSchemes.path()
            << ".  Copy from base" << endl;

        IOdictionary fvSchemesBase
        (
            IOobject
            (
                "fvSchemes",
                baseMesh().time().system(),
                baseMesh().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        fvSchemes = fvSchemesBase;
        fvSchemes.regIOobject::write();
    }

    IOdictionary fvSolution
    (
        IOobject
        (
            "fvSolution",
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    if (!fvSolution.headerOk())
    {
        Info << "Cannot read " << fvSolution.path()
            << ".  Copy from base" << endl;

        IOdictionary fvSolutionBase
        (
            IOobject
            (
                "fvSolution",
                baseMesh().time().system(),
                baseMesh().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        fvSolution = fvSolutionBase;
        fvSolution.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::newFvMeshSubset::newFvMeshSubset
(
    const IOobject& io,
    const fvMesh& baseMesh
)
:
    regIOobject(io),
    baseMesh_(baseMesh),
    newFvMeshSubsetPtr_(NULL),
    pointMeshSubsetPtr_(NULL),
    pointMap_(0),
    faceMap_(0),
    cellMap_(0),
    patchMap_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newFvMeshSubset::~newFvMeshSubset()
{
    deleteDemandDrivenData(pointMeshSubsetPtr_);
    deleteDemandDrivenData(newFvMeshSubsetPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::newFvMeshSubset::setCellSubset
(
    const labelHashSet& globalCellMap,
    const label patchID
)
{
    // Initial check on patches before doing anything time consuming.
    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();
    const cellList& oldCells = baseMesh().cells();
    const faceList& oldFaces = baseMesh().allFaces();
    const pointField& oldPoints = baseMesh().points();
    const labelList& oldOwner = baseMesh().faceOwner();
    const labelList& oldNeighbour = baseMesh().faceNeighbour();

    label wantedPatchID = patchID;

    if (wantedPatchID == -1)
    {
        // No explicit patch specified. Put in oldInternalFaces patch.
        // Check if patch with this name already exists.
        wantedPatchID = oldPatches.findPatchID("oldInternalFaces");
    }
    else if (wantedPatchID < 0 || wantedPatchID >= oldPatches.size())
    {
        FatalErrorIn
        (
            "newFvMeshSubset::setCellSubset(const labelHashSet&"
            ", const label patchID)"
        )   << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << oldPatches.size()-1
            << abort(FatalError);
    }


    cellMap_ = globalCellMap.toc();

    // Sort the cell map in the ascending order
    sort(cellMap_);

    label nCellsInSet = cellMap_.size();

    // Mark all used faces

    Map<label> facesToSubset(primitiveMesh::facesPerCell_*nCellsInSet);

    forAll (cellMap_, cellI)
    {
        // Mark all faces from the cell
        const labelList& curFaces = oldCells[cellMap_[cellI]];

        forAll (curFaces, faceI)
        {
            if (!facesToSubset.found(curFaces[faceI]))
            {
                facesToSubset.insert(curFaces[faceI], 1);
            }
            else
            {
                facesToSubset[curFaces[faceI]]++;
            }
        }
    }

    // Mark all used points and make a global-to-local face map
    Map<label> globalFaceMap(facesToSubset.size());

    // Make a global-to-local point map
    Map<label> globalPointMap
    (
        primitiveMesh::pointsPerFace_*facesToSubset.size()
    );

    // This is done in two goes, so that the boundary faces are last
    // in the list.  Because of this, I need to create the face map
    // along the way rather than just grab the table of contents.
    labelList facesToc = facesToSubset.toc();
    sort(facesToc);
    faceMap_.setSize(facesToc.size());

    // 1. Get all faces that will be internal to the submesh.
    forAll (facesToc, faceI)
    {
        if (facesToSubset[facesToc[faceI]] == 2)
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[faceI];
            globalFaceMap.insert(facesToc[faceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[faceI]], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    label nInternalFaces = globalFaceMap.size();


    // Where to insert old internal faces.
    label oldPatchStart = labelMax;
    if (wantedPatchID != -1)
    {
        oldPatchStart = oldPatches[wantedPatchID].start();
    }


    label faceI = 0;

    // 2. Boundary faces up to where we want to insert old internal faces
    for (; faceI< facesToc.size(); faceI++)
    {
        if (facesToc[faceI] >= oldPatchStart)
        {
            break;
        }
        if
        (
            !baseMesh().isInternalFace(facesToc[faceI])
         && facesToSubset[facesToc[faceI]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[faceI];
            globalFaceMap.insert(facesToc[faceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[faceI]], globalPointMap);
        }
    }

    // 3. old internal faces
    forAll(facesToc, intFaceI)
    {
        if
        (
            baseMesh().isInternalFace(facesToc[intFaceI])
         && facesToSubset[facesToc[intFaceI]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[intFaceI];
            globalFaceMap.insert(facesToc[intFaceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[intFaceI]], globalPointMap);
        }
    }

    // 4. Remaining boundary faces
    for (; faceI< facesToc.size(); faceI++)
    {
        if
        (
            !baseMesh().isInternalFace(facesToc[faceI])
         && facesToSubset[facesToc[faceI]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[faceI];
            globalFaceMap.insert(facesToc[faceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[faceI]], globalPointMap);
        }
    }



    // Grab the points map
    pointMap_ = globalPointMap.toc();
    sort(pointMap_);

    forAll (pointMap_, pointI)
    {
        globalPointMap[pointMap_[pointI]] = pointI;
    }

    if (fvMesh::debug)
    {
        Pout<< "Number of cells in new mesh: " << nCellsInSet << nl
            << "Number of faces in new mesh: " << globalFaceMap.size() << nl
            << "Number of points in new mesh: " << globalPointMap.size()
            << endl;
    }

    // Make a new mesh
    pointField newPoints(globalPointMap.size());

    label nNewPoints = 0;

    forAll (pointMap_, pointI)
    {
        newPoints[nNewPoints] = oldPoints[pointMap_[pointI]];
        nNewPoints++;
    }

    faceList newFaces(globalFaceMap.size());

    label nNewFaces = 0;

    // Make internal faces
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const face& oldF = oldFaces[faceMap_[faceI]];

        face newF(oldF.size());

        forAll (newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }

    // Make boundary faces

    label nbSize = oldPatches.size();
    label oldInternalPatchID  = -1;

    if (wantedPatchID == -1)
    {
        // Create 'oldInternalFaces' patch at the end
        // and put all exposed internal faces in there.
        oldInternalPatchID = nbSize;
        nbSize++;

    }
    else
    {
        oldInternalPatchID = wantedPatchID;
    }


    // Grad size and start of each patch on the fly.  Because of the
    // structure of the underlying mesh, the patches will appear in the
    // ascending order
    labelList boundaryPatchSizes(nbSize, 0);

    // Assign boundary faces. Visited in order of faceMap_.
    for (label faceI = nInternalFaces; faceI < faceMap_.size(); faceI++)
    {
        label oldFaceI = faceMap_[faceI];

        face oldF = oldFaces[oldFaceI];

        // Turn the faces as necessary to point outwards
        if (baseMesh().isInternalFace(oldFaceI))
        {
            // Internal face. Possibly turned the wrong way round
            if
            (
                !globalCellMap.found(oldOwner[oldFaceI])
             && globalCellMap.found(oldNeighbour[oldFaceI])
            )
            {
                oldF = oldFaces[oldFaceI].reverseFace();
            }

            // Update count for patch
            boundaryPatchSizes[oldInternalPatchID]++;
        }
        else
        {
            // Boundary face. Increment the appropriate patch
            label patchOfFace = oldPatches.whichPatch(oldFaceI);

            // Update count for patch
            boundaryPatchSizes[patchOfFace]++;
        }

        face newF(oldF.size());

        forAll (newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }



    // Create cells
    cellList newCells(nCellsInSet);

    label nNewCells = 0;

    forAll (cellMap_, cellI)
    {
        const labelList& oldC = oldCells[cellMap_[cellI]];

        labelList newC(oldC.size());

        forAll (newC, i)
        {
            newC[i] = globalFaceMap[oldC[i]];
        }

        newCells[nNewCells] = cell(newC);
        nNewCells++;
    }


    // Make a new mesh
    newFvMeshSubsetPtr_ = new fvMesh
    (
        IOobject
        (
            this->name() + "Subset",
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newCells)
    );

    // Clear point mesh
    deleteDemandDrivenData(pointMeshSubsetPtr_);


    // Add old patches
    List<polyPatch*> newBoundary(nbSize);
    patchMap_.setSize(nbSize);
    label nNewPatches = 0;
    label patchStart = nInternalFaces;

    forAll(oldPatches, patchI)
    {
        if (boundaryPatchSizes[patchI] > 0)
        {
            // Patch still exists. Add it
            newBoundary[nNewPatches] = oldPatches[patchI].clone
            (
                newFvMeshSubsetPtr_->boundaryMesh(),
                nNewPatches,
                boundaryPatchSizes[patchI],
                patchStart
            ).ptr();

            patchStart += boundaryPatchSizes[patchI];
            patchMap_[nNewPatches] = patchI;
            nNewPatches++;
        }
    }

    if (wantedPatchID == -1)
    {
        // Newly created patch so is at end. Check if any faces in it.
        if (boundaryPatchSizes[oldInternalPatchID] > 0)
        {
            newBoundary[nNewPatches] = new polyPatch
            (
                "oldInternalFaces",
                boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                newFvMeshSubsetPtr_->boundaryMesh()
            );

            // The index for the first patch is -1 as it originates from
            // the internal faces
            patchMap_[nNewPatches] = -1;
            nNewPatches++;
        }
    }

    // Reset the patch lists
    newBoundary.setSize(nNewPatches);
    patchMap_.setSize(nNewPatches);

    // Add the fvPatches
    newFvMeshSubsetPtr_->addFvPatches(newBoundary);

    // Subset and add any zones
    subsetZones();
}


void Foam::newFvMeshSubset::setLargeCellSubset
(
    const labelList& region,
    const label currentRegion,
    const label patchID,
    const bool syncPar
)
{
    const cellList& oldCells = baseMesh().cells();
    const faceList& oldFaces = baseMesh().faces();
    const pointField& oldPoints = baseMesh().points();
    const labelList& oldOwner = baseMesh().faceOwner();
    const labelList& oldNeighbour = baseMesh().faceNeighbour();
    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();
    const label oldNInternalFaces = baseMesh().nInternalFaces();

    // Initial checks

    if (region.size() != oldCells.size())
    {
        FatalErrorIn
        (
            "newFvMeshSubset::setCellSubset(const labelList&"
            ", const label, const label, const bool)"
        )   << "Size of region " << region.size()
            << " is not equal to number of cells in mesh " << oldCells.size()
            << abort(FatalError);
    }


    label wantedPatchID = patchID;

    if (wantedPatchID == -1)
    {
        // No explicit patch specified. Put in oldInternalFaces patch.
        // Check if patch with this name already exists.
        wantedPatchID = oldPatches.findPatchID("oldInternalFaces");
    }
    else if (wantedPatchID < 0 || wantedPatchID >= oldPatches.size())
    {
        FatalErrorIn
        (
            "newFvMeshSubset::setCellSubset(const labelList&"
            ", const label, const label, const bool)"
        )   << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << oldPatches.size()-1
            << abort(FatalError);
    }


    // Get the cells for the current region.
    cellMap_.setSize(oldCells.size());
    label nCellsInSet = 0;

    forAll (region, oldCellI)
    {
        if (region[oldCellI] == currentRegion)
        {
            cellMap_[nCellsInSet++] = oldCellI;
        }
    }
    cellMap_.setSize(nCellsInSet);


    // Mark all used faces. Count number of cells using them
    // 0: face not used anymore
    // 1: face used by one cell, face becomes/stays boundary face
    // 2: face still used and remains internal face
    // 3: face coupled and used by one cell only (so should become normal,
    //    non-coupled patch face)
    //
    // Note that this is not really nessecary - but means we can size things
    // correctly. Also makes handling coupled faces much easier.

    labelList nCellsUsingFace(oldFaces.size(), 0);

    label nFacesInSet = 0;
    forAll (oldFaces, oldFaceI)
    {
        bool faceUsed = false;

        if (region[oldOwner[oldFaceI]] == currentRegion)
        {
            nCellsUsingFace[oldFaceI]++;
            faceUsed = true;
        }

        if
        (
            baseMesh().isInternalFace(oldFaceI)
         && (region[oldNeighbour[oldFaceI]] == currentRegion)
        )
        {
            nCellsUsingFace[oldFaceI]++;
            faceUsed = true;
        }

        if (faceUsed)
        {
            nFacesInSet++;
        }
    }
    faceMap_.setSize(nFacesInSet);

    // Handle coupled faces. Modifies patch faces to be uncoupled to 3.
    doCoupledPatches(syncPar, nCellsUsingFace);


    // See which patch to use for exposed internal faces.
    label oldInternalPatchID = 0;

    // Insert faces before which patch
    label nextPatchID = oldPatches.size();

    // old to new patches
    labelList globalPatchMap(oldPatches.size());

    // New patch size
    label nbSize = oldPatches.size();

    if (wantedPatchID == -1)
    {
        // Create 'oldInternalFaces' patch at the end (or before
        // processorPatches)
        // and put all exposed internal faces in there.

        forAll(oldPatches, patchI)
        {
            if (isA<processorPolyPatch>(oldPatches[patchI]))
            {
                nextPatchID = patchI;
                break;
            }
            oldInternalPatchID++;
        }

        nbSize++;

        // adapt old to new patches for inserted patch
        for (label oldPatchI = 0; oldPatchI < nextPatchID; oldPatchI++)
        {
            globalPatchMap[oldPatchI] = oldPatchI;
        }
        for
        (
            label oldPatchI = nextPatchID;
            oldPatchI < oldPatches.size();
            oldPatchI++
        )
        {
            globalPatchMap[oldPatchI] = oldPatchI+1;
        }
    }
    else
    {
        oldInternalPatchID = wantedPatchID;
        nextPatchID = wantedPatchID+1;

        // old to new patches
        globalPatchMap = identity(oldPatches.size());
    }

    labelList boundaryPatchSizes(nbSize, 0);


    // Make a global-to-local point map
    labelList globalPointMap(oldPoints.size(), -1);

    labelList globalFaceMap(oldFaces.size(), -1);
    label faceI = 0;

    // 1. Pick up all preserved internal faces.
    for (label oldFaceI = 0; oldFaceI < oldNInternalFaces; oldFaceI++)
    {
        if (nCellsUsingFace[oldFaceI] == 2)
        {
            globalFaceMap[oldFaceI] = faceI;
            faceMap_[faceI++] = oldFaceI;

            // Mark all points from the face
            markPoints(oldFaces[oldFaceI], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    label nInternalFaces = faceI;

    // 2. Boundary faces up to where we want to insert old internal faces
    for
    (
        label oldPatchI = 0;
        oldPatchI < oldPatches.size()
     && oldPatchI < nextPatchID;
        oldPatchI++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchI];

        label oldFaceI = oldPatch.start();

        forAll (oldPatch, i)
        {
            if (nCellsUsingFace[oldFaceI] == 1)
            {
                // Boundary face is kept.

                // Mark face and increment number of points in set
                globalFaceMap[oldFaceI] = faceI;
                faceMap_[faceI++] = oldFaceI;

                // Mark all points from the face
                markPoints(oldFaces[oldFaceI], globalPointMap);

                // Increment number of patch faces
                boundaryPatchSizes[globalPatchMap[oldPatchI]]++;
            }
            oldFaceI++;
        }
    }

    // 3a. old internal faces that have become exposed.
    for (label oldFaceI = 0; oldFaceI < oldNInternalFaces; oldFaceI++)
    {
        if (nCellsUsingFace[oldFaceI] == 1)
        {
            globalFaceMap[oldFaceI] = faceI;
            faceMap_[faceI++] = oldFaceI;

            // Mark all points from the face
            markPoints(oldFaces[oldFaceI], globalPointMap);

            // Increment number of patch faces
            boundaryPatchSizes[oldInternalPatchID]++;
        }
    }

    // 3b. coupled patch faces that have become uncoupled.
    for
    (
        label oldFaceI = oldNInternalFaces;
        oldFaceI < oldFaces.size();
        oldFaceI++
    )
    {
        if (nCellsUsingFace[oldFaceI] == 3)
        {
            globalFaceMap[oldFaceI] = faceI;
            faceMap_[faceI++] = oldFaceI;

            // Mark all points from the face
            markPoints(oldFaces[oldFaceI], globalPointMap);

            // Increment number of patch faces
            boundaryPatchSizes[oldInternalPatchID]++;
        }
    }

    // 4. Remaining boundary faces
    for
    (
        label oldPatchI = nextPatchID;
        oldPatchI < oldPatches.size();
        oldPatchI++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchI];

        label oldFaceI = oldPatch.start();

        forAll (oldPatch, i)
        {
            if (nCellsUsingFace[oldFaceI] == 1)
            {
                // Boundary face is kept.

                // Mark face and increment number of points in set
                globalFaceMap[oldFaceI] = faceI;
                faceMap_[faceI++] = oldFaceI;

                // Mark all points from the face
                markPoints(oldFaces[oldFaceI], globalPointMap);

                // Increment number of patch faces
                boundaryPatchSizes[globalPatchMap[oldPatchI]]++;
            }
            oldFaceI++;
        }
    }

    if (faceI != nFacesInSet)
    {
        FatalErrorIn
        (
            "newFvMeshSubset::setCellSubset(const labelList&"
            ", const label, const label, const bool)"
        )   << "Problem" << abort(FatalError);
    }


    // Grab the points map
    label nPointsInSet = 0;

    forAll(globalPointMap, pointI)
    {
        if (globalPointMap[pointI] != -1)
        {
            nPointsInSet++;
        }
    }
    pointMap_.setSize(nPointsInSet);

    nPointsInSet = 0;

    forAll(globalPointMap, pointI)
    {
        if (globalPointMap[pointI] != -1)
        {
            pointMap_[nPointsInSet] = pointI;
            globalPointMap[pointI] = nPointsInSet;
            nPointsInSet++;
        }
    }

    Pout<< "Number of cells in new mesh : " << cellMap_.size() << endl;
    Pout<< "Number of faces in new mesh : " << faceMap_.size() << endl;
    Pout<< "Number of points in new mesh: " << pointMap_.size() << endl;

    // Make a new mesh
    pointField newPoints(pointMap_.size());

    label nNewPoints = 0;

    forAll (pointMap_, pointI)
    {
        newPoints[nNewPoints] = oldPoints[pointMap_[pointI]];
        nNewPoints++;
    }

    faceList newFaces(faceMap_.size());

    label nNewFaces = 0;

    // Make internal faces
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const face& oldF = oldFaces[faceMap_[faceI]];

        face newF(oldF.size());

        forAll (newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }


    // Make boundary faces. (different from internal since might need to be
    // flipped)
    for (label faceI = nInternalFaces; faceI < faceMap_.size(); faceI++)
    {
        label oldFaceI = faceMap_[faceI];

        face oldF = oldFaces[oldFaceI];

        // Turn the faces as necessary to point outwards
        if (baseMesh().isInternalFace(oldFaceI))
        {
            // Was internal face. Possibly turned the wrong way round
            if
            (
                region[oldOwner[oldFaceI]] != currentRegion
             && region[oldNeighbour[oldFaceI]] == currentRegion
            )
            {
                oldF = oldFaces[oldFaceI].reverseFace();
            }
        }

        // Relabel vertices of the (possibly turned) face.
        face newF(oldF.size());

        forAll (newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }



    // Create cells
    cellList newCells(nCellsInSet);

    label nNewCells = 0;

    forAll (cellMap_, cellI)
    {
        const labelList& oldC = oldCells[cellMap_[cellI]];

        labelList newC(oldC.size());

        forAll (newC, i)
        {
            newC[i] = globalFaceMap[oldC[i]];
        }

        newCells[nNewCells] = cell(newC);
        nNewCells++;
    }


    // Make a new mesh
    makeFvDictionaries();

    newFvMeshSubsetPtr_ = new fvMesh
    (
        IOobject
        (
            this->name() + "Subset",
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newCells),
        syncPar           // parallel synchronisation
    );

    // Clear point mesh
    deleteDemandDrivenData(pointMeshSubsetPtr_);

    // Add old patches
    List<polyPatch*> newBoundary(nbSize);
    patchMap_.setSize(nbSize);
    label nNewPatches = 0;
    label patchStart = nInternalFaces;


    // For parallel: only remove patch if none of the processors has it.
    // This only gets done for patches before the one being inserted
    // (so patches < nextPatchID)

    // Get sum of patch sizes. Zero if patch can be deleted.
    labelList globalPatchSizes(boundaryPatchSizes);
    globalPatchSizes.setSize(nextPatchID);

    if (syncPar && Pstream::parRun())
    {
        // Get patch names (up to nextPatchID)
        List<wordList> patchNames(Pstream::nProcs());
        patchNames[Pstream::myProcNo()] = oldPatches.names();
        patchNames[Pstream::myProcNo()].setSize(nextPatchID);
        Pstream::gatherList(patchNames);
        Pstream::scatterList(patchNames);

        // Get patch sizes (up to nextPatchID).
        // Note that up to nextPatchID the globalPatchMap is an identity so
        // no need to index through that.
        Pstream::listCombineGather(globalPatchSizes, plusEqOp<label>());
        Pstream::listCombineScatter(globalPatchSizes);

        // Now all processors have all the patchnames.
        // Decide: if all processors have the same patch names and size is zero
        // everywhere remove the patch.
        bool samePatches = true;

        for (label procI = 1; procI < patchNames.size(); procI++)
        {
            if (patchNames[procI] != patchNames[0])
            {
                samePatches = false;
                break;
            }
        }

        if (!samePatches)
        {
            // Patchnames not sync on all processors so disable removal of
            // zero sized patches.
            globalPatchSizes = labelMax;
        }
    }


    // Old patches

    for
    (
        label oldPatchI = 0;
        oldPatchI < oldPatches.size()
     && oldPatchI < nextPatchID;
        oldPatchI++
    )
    {
        label newSize = boundaryPatchSizes[globalPatchMap[oldPatchI]];

#if FOAMEXTEND > 40
        if (isA<oversetPolyPatch>(oldPatches[oldPatchI]))
        {
            newBoundary[nNewPatches] = new polyPatch
            (
                oldPatches[oldPatchI].name(),
                newSize, //boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                newFvMeshSubsetPtr_->boundaryMesh()
            );
        }
        else
#endif
        {
            // Clone (even if 0 size)
            newBoundary[nNewPatches] = oldPatches[oldPatchI].clone
            (
                newFvMeshSubsetPtr_->boundaryMesh(),
                nNewPatches,
                newSize,
                patchStart
            ).ptr();
        }

        patchStart += newSize;
        patchMap_[nNewPatches] = oldPatchI;    // compact patchMap
        nNewPatches++;
    }

    // Inserted patch

    if (wantedPatchID == -1)
    {
        label oldInternalSize = boundaryPatchSizes[oldInternalPatchID];

        if (syncPar)
        {
            reduce(oldInternalSize, sumOp<label>());
        }

        // Newly created patch so is at end. Check if any faces in it.
        if (oldInternalSize > 0)
        {
            newBoundary[nNewPatches] = new polyPatch
            (
                "oldInternalFaces",
                boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                newFvMeshSubsetPtr_->boundaryMesh()
            );

            Pout<< "    oldInternalFaces : "
                << boundaryPatchSizes[oldInternalPatchID] << endl;

            // The index for the first patch is -1 as it originates from
            // the internal faces
            patchStart += boundaryPatchSizes[oldInternalPatchID];
            patchMap_[nNewPatches] = -1;
            nNewPatches++;
        }
    }

    // Old patches

    for
    (
        label oldPatchI = nextPatchID;
        oldPatchI < oldPatches.size();
        oldPatchI++
    )
    {
        label newSize = boundaryPatchSizes[globalPatchMap[oldPatchI]];

        // Patch still exists. Add it
        newBoundary[nNewPatches] = oldPatches[oldPatchI].clone
        (
            newFvMeshSubsetPtr_->boundaryMesh(),
            nNewPatches,
            newSize,
            patchStart
        ).ptr();

        //Pout<< "    " << oldPatches[oldPatchI].name() << " : "
        //    << newSize << endl;

        patchStart += newSize;
        patchMap_[nNewPatches] = oldPatchI;    // compact patchMap
        nNewPatches++;
    }

    // Reset the patch lists
    newBoundary.setSize(nNewPatches);
    patchMap_.setSize(nNewPatches);

    // Add the fvPatches
    newFvMeshSubsetPtr_->addFvPatches(newBoundary, syncPar);

    // Subset and add any zones
    subsetZones();
}


void Foam::newFvMeshSubset::setLargeCellSubset
(
    const labelHashSet& globalCellMap,
    const label patchID,
    const bool syncPar
)
{
    labelList region(baseMesh().nCells(), 0);

    forAllConstIter (labelHashSet, globalCellMap, iter)
    {
        region[iter.key()] = 1;
    }
    setLargeCellSubset(region, 1, patchID, syncPar);
}


const fvMesh& Foam::newFvMeshSubset::subMesh() const
{
    checkCellSubset();

    return *newFvMeshSubsetPtr_;
}


fvMesh& Foam::newFvMeshSubset::subMesh()
{
    checkCellSubset();

    return *newFvMeshSubsetPtr_;
}


const pointMesh& Foam::newFvMeshSubset::subPointMesh() const
{
    if (!pointMeshSubsetPtr_)
    {
        pointMeshSubsetPtr_ = new pointMesh(subMesh());

//         if (Pstream::parRun())
//         {
//             pointMeshSubsetPtr_ = new pointMesh(subMesh(), true);
//         }
//         else
//         {
//             pointMeshSubsetPtr_ = new pointMesh(subMesh());
//         }
    }

    return *pointMeshSubsetPtr_;
}


pointMesh& Foam::newFvMeshSubset::subPointMesh()
{
    if (!pointMeshSubsetPtr_)
    {
        pointMeshSubsetPtr_ = new pointMesh(subMesh(), true);
    }

    return *pointMeshSubsetPtr_;
}


const labelList& Foam::newFvMeshSubset::pointMap() const
{
    checkCellSubset();

    return pointMap_;
}


const labelList& Foam::newFvMeshSubset::faceMap() const
{
    checkCellSubset();

    return faceMap_;
}


const labelList& Foam::newFvMeshSubset::cellMap() const
{
    checkCellSubset();

    return cellMap_;
}


const labelList& Foam::newFvMeshSubset::patchMap() const
{
    checkCellSubset();

    return patchMap_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
