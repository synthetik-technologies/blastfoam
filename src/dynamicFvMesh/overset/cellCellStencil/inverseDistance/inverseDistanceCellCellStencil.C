/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "inverseDistanceCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
#include "Time.H"
#include "fvMeshSubset.H"

#include "globalIndex.H"
#include "oversetFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "treeBoundBoxList.H"
#include "waveMethod.H"

#include "regionSplit.H"
#include "dynamicOversetBlastFvMesh.H"
#include "OSspecific.H"
//#include "minData.H"
//#include "FaceCellWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(inverseDistance, 0);
    addToRunTimeSelectionTable(cellCellStencil, inverseDistance, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::cellCellStencils::inverseDistance::index
(
    const labelVector& nDivs,
    const labelVector& ijk
)
{
    return (ijk[0]*nDivs[1] + ijk[1])*nDivs[2] + ijk[2];
}


Foam::labelVector Foam::cellCellStencils::inverseDistance::index3
(
    const labelVector& nDivs,
    const label boxI
)
{
    label ij = boxI/nDivs[2];
    label k = boxI-ij*nDivs[2];
    label i = ij/nDivs[1];
    label j = ij-i*nDivs[1];

    return labelVector(i, j, k);
}


Foam::labelVector Foam::cellCellStencils::inverseDistance::index3
(
    const boundBox& bb,
    const labelVector& nDivs,
    const point& pt
)
{
    const vector d(bb.span());
    const point relPt(pt-bb.min());

    return labelVector
    (
        floor(relPt[0]/d[0]*nDivs[0]),
        floor(relPt[1]/d[1]*nDivs[1]),
        floor(relPt[2]/d[2]*nDivs[2])
    );
}


Foam::point Foam::cellCellStencils::inverseDistance::position
(
    const boundBox& bb,
    const labelVector& nDivs,
    const label boxI
)
{
    // Return midpoint of box indicated by boxI
    labelVector ids(index3(nDivs, boxI));

    const vector d(bb.span());
    const vector sz(d[0]/nDivs[0], d[1]/nDivs[1], d[2]/nDivs[2]);

    return bb.min()+0.5*sz+vector(sz[0]*ids[0], sz[1]*ids[1], sz[2]*ids[2]);
}


void Foam::cellCellStencils::inverseDistance::fill
(
    PackedList<2>& elems,
    const boundBox& bb,
    const labelVector& nDivs,
    const boundBox& subBb,
    const unsigned int val
)
{
    labelVector minIds(index3(bb, nDivs, subBb.min()));
    labelVector maxIds(index3(bb, nDivs, subBb.max()));

    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        if (maxIds[cmpt] < 0 || minIds[cmpt] > nDivs[cmpt])
        {
            return;
        }
    }

    labelVector maxIndex(labelVector(nDivs[0]-1, nDivs[1]-1, nDivs[2]-1));
    minIds = max(labelVector::zero, minIds);
    maxIds = min(maxIndex, maxIds);

    for (label i = minIds[0]; i <= maxIds[0]; i++)
    {
        for (label j = minIds[1]; j <= maxIds[1]; j++)
        {
            for (label k = minIds[2]; k <= maxIds[2]; k++)
            {
                label i1 = index(nDivs, labelVector(i, j, k));
                elems[i1] = val;
            }
        }
    }
}


void Foam::cellCellStencils::inverseDistance::markBoundaries
(
    const fvMesh& mesh,
    const vector& smallVec,

    const boundBox& bb,
    const labelVector& nDivs,
    PackedList<2>& patchTypes,

    const labelList& cellMap,
    labelList& patchCellTypes
)
{
    // Mark all voxels that overlap the bounding box of any patch

    const fvBoundaryMesh& pbm = mesh.boundary();

    patchTypes = patchCellType::OTHER;

    // Mark wall boundaries
    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];
        const labelList& fc = fvp.faceCells();

        if (!fvPatch::constraintType(fvp.type()))
        {
            //Info<< "Marking cells on proper patch " << fvp.name()
            //    << " with type " << fvp.type() << endl;
            const polyPatch& pp = fvp.patch();
            forAll(pp, i)
            {
                // Mark in overall patch types
                patchCellTypes[cellMap[fc[i]]] = patchCellType::PATCH;

                // Mark in voxel mesh
                boundBox faceBb(pp.points(), pp[i], false);
                faceBb.min() -= smallVec;
                faceBb.max() += smallVec;

                if (bb.overlaps(faceBb))
                {
                    fill(patchTypes, bb, nDivs, faceBb, patchCellType::PATCH);
                }
            }
        }
    }

    // Override with overset boundaries
    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];
        const labelList& fc = fvp.faceCells();

        if (isA<oversetFvPatch>(fvp))
        {
            //Info<< "Marking cells on overset patch " << fvp.name() << endl;
            const polyPatch& pp = fvp.patch();
            forAll(pp, i)
            {
                // Mark in overall patch types
                patchCellTypes[cellMap[fc[i]]] = patchCellType::OVERSET;

                // Mark in voxel mesh
                boundBox faceBb(pp.points(), pp[i], false);
                faceBb.min() -= smallVec;
                faceBb.max() += smallVec;

                if (bb.overlaps(faceBb))
                {
                    fill(patchTypes, bb, nDivs, faceBb, patchCellType::OVERSET);
                }
            }
        }
    }
}


Foam::treeBoundBox Foam::cellCellStencils::inverseDistance::cellBb
(
    const primitiveMesh& mesh,
    const label celli
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    treeBoundBox bb
    (
        vector(GREAT, GREAT, GREAT),
        vector(-GREAT, -GREAT, -GREAT)
    );

    const cell& cFaces = cells[celli];

    forAll(cFaces, cFacei)
    {
        const face& f = faces[cFaces[cFacei]];

        forAll(f, fp)
        {
            const point& p = points[f[fp]];

            bb.min() = min(bb.min(), p);
            bb.max() = max(bb.max(), p);
        }
    }
    return bb;
}


bool Foam::cellCellStencils::inverseDistance::overlaps
(
    const boundBox& bb,
    const labelVector& nDivs,
    const PackedList<2>& vals,
    const treeBoundBox& subBb,
    const unsigned int val
)
{
    // Checks if subBb overlaps any voxel set to val

    labelVector minIds(index3(bb, nDivs, subBb.min()));
    labelVector maxIds(index3(bb, nDivs, subBb.max()));

    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        if (maxIds[cmpt] < 0 || minIds[cmpt] > nDivs[cmpt])
        {
            return false;
        }
    }

    labelVector maxIndex(labelVector(nDivs[0]-1, nDivs[1]-1, nDivs[2]-1));
    minIds = max(labelVector::zero, minIds);
    maxIds = min(maxIndex, maxIds);

    for (label i = minIds[0]; i <= maxIds[0]; i++)
    {
        for (label j = minIds[1]; j <= maxIds[1]; j++)
        {
            for (label k = minIds[2]; k <= maxIds[2]; k++)
            {
                label i1 = index(nDivs, labelVector(i, j, k));
                if (vals[i1] == patchCellType::PATCH)
                {
                    return true;
                }
            }
        }
    }
    return false;
}


void Foam::cellCellStencils::inverseDistance::markPatchesAsHoles
(
    PstreamBuffers& pBufs,

    const PtrList<fvMeshSubset>& meshParts,

    const List<treeBoundBoxList>& patchBb,
    const List<labelVector>& patchDivisions,
    const PtrList<PackedList<2>>& patchParts,

    const label srcI,
    const label tgtI,
    labelList& allCellTypes
) const
{
    const treeBoundBoxList& srcPatchBbs = patchBb[srcI];
    const treeBoundBoxList& tgtPatchBbs = patchBb[tgtI];
    const labelList& tgtCellMap = meshParts[tgtI].cellMap();

    // 1. do processor-local src-tgt patch overlap
    {
        const treeBoundBox& srcPatchBb = srcPatchBbs[Pstream::myProcNo()];
        const treeBoundBox& tgtPatchBb = tgtPatchBbs[Pstream::myProcNo()];

        if (srcPatchBb.overlaps(tgtPatchBb))
        {
            const PackedList<2>& srcPatchTypes = patchParts[srcI];
            const labelVector& zoneDivs = patchDivisions[srcI];

            forAll(tgtCellMap, tgtCelli)
            {
                label celli = tgtCellMap[tgtCelli];
                treeBoundBox cBb(cellBb(mesh_, celli));
                cBb.min() -= smallVec_;
                cBb.max() += smallVec_;

                if
                (
                    overlaps
                    (
                        srcPatchBb,
                        zoneDivs,
                        srcPatchTypes,
                        cBb,
                        patchCellType::PATCH
                    )
                )
                {
                    allCellTypes[celli] = HOLE;
                }
            }
        }
    }


    // 2. Send over srcMesh bits that overlap tgt and do calculation
    pBufs.clear();
//     for (const int procI : Pstream::allProcs())
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            const treeBoundBox& srcPatchBb = srcPatchBbs[Pstream::myProcNo()];
            const treeBoundBox& tgtPatchBb = tgtPatchBbs[procI];

            if (srcPatchBb.overlaps(tgtPatchBb))
            {
                // Send over complete patch voxel map. Tbd: could
                // subset
                UOPstream os(procI, pBufs);
                os << srcPatchBb << patchDivisions[srcI] << patchParts[srcI];
            }
        }
    }
    pBufs.finishedSends();
//     for (const int procI : Pstream::allProcs())
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            //const treeBoundBox& srcBb = srcBbs[procI];
            const treeBoundBox& srcPatchBb = srcPatchBbs[procI];
            const treeBoundBox& tgtPatchBb = tgtPatchBbs[Pstream::myProcNo()];

            if (srcPatchBb.overlaps(tgtPatchBb))
            {
                UIPstream is(procI, pBufs);
                {
                    treeBoundBox receivedBb(is);
                    if (srcPatchBb != receivedBb)
                    {
                        FatalErrorInFunction
                            << "proc:" << procI
                            << " srcPatchBb:" << srcPatchBb
                            << " receivedBb:" << receivedBb
                            << exit(FatalError);
                    }
                }
                const labelVector zoneDivs(is);
                const PackedList<2> srcPatchTypes(is);

                forAll(tgtCellMap, tgtCelli)
                {
                    label celli = tgtCellMap[tgtCelli];
                    treeBoundBox cBb(cellBb(mesh_, celli));
                    cBb.min() -= smallVec_;
                    cBb.max() += smallVec_;
                    if
                    (
                        overlaps
                        (
                            srcPatchBb,
                            zoneDivs,
                            srcPatchTypes,
                            cBb,
                            patchCellType::PATCH
                        )
                    )
                    {
                        allCellTypes[celli] = HOLE;
                    }
                }
            }
        }
    }
}


bool Foam::cellCellStencils::inverseDistance::betterDonor
(
    const label destMesh,
    const label currentDonorMesh,
    const label newDonorMesh
) const
{
    // This determines for multiple overlapping meshes which one provides
    // the best donors. Is very basic and only looks at indices of meshes:
    // - 'nearest' mesh index wins, i.e. on mesh 0 it preferentially uses donors
    //   from mesh 1 over mesh 2 (if applicable)
    // - if same 'distance' the highest mesh wins. So on mesh 1 it
    //   preferentially uses donors from mesh 2 over mesh 0. This particular
    //   rule helps to avoid some interpolation loops where mesh 1 uses donors
    //   from mesh 0 (usually the background) but mesh 0 then uses
    //   donors from 1.

    if (currentDonorMesh == -1)
    {
        return true;
    }
    else
    {
        const label currentDist = mag(currentDonorMesh-destMesh);
        const label newDist = mag(newDonorMesh-destMesh);

        if (newDist < currentDist)
        {
            return true;
        }
        else if (newDist == currentDist && newDonorMesh > currentDonorMesh)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


void Foam::cellCellStencils::inverseDistance::markDonors
(
    const globalIndex& globalCells,
    PstreamBuffers& pBufs,
    const PtrList<fvMeshSubset>& meshParts,
    const List<treeBoundBoxList>& meshBb,

    const labelList& allCellTypes,

    const label srcI,
    const label tgtI,
    labelListList& allStencil,
    labelList& allDonor
) const
{
    const treeBoundBoxList& srcBbs = meshBb[srcI];
    const treeBoundBoxList& tgtBbs = meshBb[tgtI];

    const fvMesh& srcMesh = meshParts[srcI].subMesh();
    const labelList& srcCellMap = meshParts[srcI].cellMap();
    const fvMesh& tgtMesh = meshParts[tgtI].subMesh();
    const pointField& tgtCc = tgtMesh.cellCentres();
    const labelList& tgtCellMap = meshParts[tgtI].cellMap();

    // 1. do processor-local src/tgt overlap
    {
        labelList tgtToSrcAddr;
        waveMethod::calculate(tgtMesh, srcMesh, tgtToSrcAddr);
        forAll(tgtCellMap, tgtCelli)
        {
            label srcCelli = tgtToSrcAddr[tgtCelli];
            if (srcCelli != -1 && allCellTypes[srcCellMap[srcCelli]] != HOLE)
            {
                label celli = tgtCellMap[tgtCelli];

                // TBD: check for multiple donors. Maybe better one? For
                //      now check 'nearer' mesh
                if (betterDonor(tgtI, allDonor[celli], srcI))
                {
                    label globalDonor =
                        globalCells.toGlobal(srcCellMap[srcCelli]);
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = globalDonor;
                    allDonor[celli] = srcI;
                }
            }
        }
    }


    // 2. Send over tgtMesh bits that overlap src and do calculation on
    //    srcMesh.


    // (remote) processors where the tgt overlaps my src
    DynamicList<label> tgtOverlapProcs(Pstream::nProcs());
    // (remote) processors where the src overlaps my tgt
    DynamicList<label> srcOverlapProcs(Pstream::nProcs());
//     for (const int procI : Pstream::allProcs())
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            if (tgtBbs[procI].overlaps(srcBbs[Pstream::myProcNo()]))
            {
                tgtOverlapProcs.append(procI);
            }
            if (srcBbs[procI].overlaps(tgtBbs[Pstream::myProcNo()]))
            {
                srcOverlapProcs.append(procI);
            }
        }
    }



    // Indices of tgtcells to send over to each processor
    List<DynamicList<label>> tgtSendCells(Pstream::nProcs());
    forAll(srcOverlapProcs, i)
    {
        label procI = srcOverlapProcs[i];
        tgtSendCells[procI].reserve(tgtMesh.nCells()/srcOverlapProcs.size());
    }


    forAll(tgtCellMap, tgtCelli)
    {
        label celli = tgtCellMap[tgtCelli];
        if (srcOverlapProcs.size())
        {
            treeBoundBox subBb(cellBb(mesh_, celli));
            subBb.min() -= smallVec_;
            subBb.max() += smallVec_;

            forAll(srcOverlapProcs, i)
            {
                label procI = srcOverlapProcs[i];
                if (subBb.overlaps(srcBbs[procI]))
                {
                    tgtSendCells[procI].append(tgtCelli);
                }
            }
        }
    }

    // Send target cell centres to overlapping processors
    pBufs.clear();

    forAll(srcOverlapProcs, i)
    {
        label procI = srcOverlapProcs[i];
        const labelList& cellIDs = tgtSendCells[procI];

        UOPstream os(procI, pBufs);
        os << UIndirectList<point>(tgtCc, cellIDs);
    }
    pBufs.finishedSends();

    // Receive bits of target processors; find; send back
    (void)srcMesh.tetBasePtIs();
    forAll(tgtOverlapProcs, i)
    {
        label procI = tgtOverlapProcs[i];

        UIPstream is(procI, pBufs);
        pointList samples(is);

        labelList donors(samples.size(), -1);
        forAll(samples, sampleI)
        {
            const point& sample = samples[sampleI];
            label srcCelli = srcMesh.findCell(sample, polyMesh::CELL_TETS);
            if (srcCelli != -1 && allCellTypes[srcCellMap[srcCelli]] != HOLE)
            {
                donors[sampleI] = globalCells.toGlobal(srcCellMap[srcCelli]);
            }
        }

        // Use same pStreamBuffers to send back.
        UOPstream os(procI, pBufs);
        os << donors;
    }
    pBufs.finishedSends();

    forAll(srcOverlapProcs, i)
    {
        label procI = srcOverlapProcs[i];
        const labelList& cellIDs = tgtSendCells[procI];

        UIPstream is(procI, pBufs);
        labelList donors(is);

        if (donors.size() != cellIDs.size())
        {
            FatalErrorInFunction<< "problem : cellIDs:" << cellIDs.size()
                << " donors:" << donors.size() << abort(FatalError);
        }

        forAll(donors, donorI)
        {
            label globalDonor = donors[donorI];

            if (globalDonor != -1)
            {
                label celli = tgtCellMap[cellIDs[donorI]];

                // TBD: check for multiple donors. Maybe better one? For
                //      now check 'nearer' mesh
                if (betterDonor(tgtI, allDonor[celli], srcI))
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = globalDonor;
                    allDonor[celli] = srcI;
                }
            }
        }
    }
}


//void Foam::cellCellStencils::inverseDistance::uncompactedRegionSplit
//(
//    const fvMesh& mesh,
//    const globalIndex& globalFaces,
//    const label nZones,
//    const labelList& zoneID,
//    const labelList& cellTypes,
//    const boolList& isBlockedFace,
//    labelList& cellRegion
//) const
//{
//    // Pass 1: locally seed 2 cells per zone (one unblocked, one blocked).
//    // This avoids excessive numbers of front
//
//    // Field on cells and faces.
//    List<minData> cellData(mesh.nCells());
//    List<minData> faceData(mesh.nFaces());
//
//    // Take over blockedFaces by seeding a negative number
//    // (so is always less than the decomposition)
//
//    forAll(isBlockedFace, facei)
//    {
//        if (isBlockedFace[facei])
//        {
//            faceData[facei] = minData(-2);
//        }
//    }
//
//
//    labelList seedFace(nZones, -1);
//
//    const labelList& owner = mesh.faceOwner();
//    const labelList& neighbour = mesh.faceNeighbour();
//
//    forAll(owner, facei)
//    {
//        label own = owner[facei];
//        if (seedFace[zoneID[own]] == -1)
//        {
//            if (cellTypes[own] != HOLE)
//            {
//                const cell& cFaces = mesh.cells()[own];
//                forAll(cFaces, i)
//                {
//                    if (!isBlockedFace[cFaces[i]])
//                    {
//                        seedFace[zoneID[own]] = cFaces[i];
//                    }
//                }
//            }
//        }
//    }
//    forAll(neighbour, facei)
//    {
//        label nei = neighbour[facei];
//        if (seedFace[zoneID[nei]] == -1)
//        {
//            if (cellTypes[nei] != HOLE)
//            {
//                const cell& cFaces = mesh.cells()[nei];
//                forAll(cFaces, i)
//                {
//                    if (!isBlockedFace[cFaces[i]])
//                    {
//                        seedFace[zoneID[nei]] = cFaces[i];
//                    }
//                }
//            }
//        }
//    }
//
//    DynamicList<label> seedFaces(nZones);
//    DynamicList<minData> seedData(seedFaces.size());
//    forAll(seedFace, zonei)
//    {
//        if (seedFace[zonei] != -1)
//        {
//            seedFaces.append(seedFace[zonei]);
//            seedData.append(minData(globalFaces.toGlobal(seedFace[zonei])));
//        }
//    }
//
//    // Propagate information inwards
//    FaceCellWave<minData> deltaCalc
//    (
//        mesh,
//        List<labelPair>(),
//        false,  // disable walking through cyclicAMI for backwards
//                // compatibility
//        seedFaces,
//        seedData,
//        faceData,
//        cellData,
//        mesh.globalData().nTotalCells()+1
//    );
//
//    // Extract
//    cellRegion.setSize(mesh.nCells());
//    forAll(cellRegion, celli)
//    {
//        if (cellData[celli].valid(deltaCalc.data()))
//        {
//            cellRegion[celli] = cellData[celli].data();
//        }
//        else
//        {
//            // Unvisited cell -> only possible if surrounded by blocked faces.
//            // If so make up region from any of the faces
//            const cell& cFaces = mesh.cells()[celli];
//            label facei = cFaces[0];
//            cellRegion[celli] = globalFaces.toGlobal(facei);
//        }
//    }
//}
//Foam::autoPtr<Foam::globalIndex>
//Foam::cellCellStencils::inverseDistance::compactedRegionSplit
//(
//    const fvMesh& mesh,
//    const globalIndex& globalRegions,
//    labelList& cellRegion
//) const
//{
//    // Now our cellRegion will have
//    // - non-local regions (i.e. originating from other processors)
//    // - non-compact locally originating regions
//    // so we'll need to compact
//
//    // 4a: count per originating processor the number of regions
//    labelList nOriginating(Pstream::nProcs(), Zero);
//    {
//        labelHashSet haveRegion(mesh.nCells()/8);
//
//        forAll(cellRegion, celli)
//        {
//            label region = cellRegion[celli];
//
//            // Count originating processor. Use isLocal as efficiency since
//            // most cells are locally originating.
//            if (globalRegions.isLocal(region))
//            {
//                if (haveRegion.insert(region))
//                {
//                    nOriginating[Pstream::myProcNo()]++;
//                }
//            }
//            else
//            {
//                label proci = globalRegions.whichProcID(region);
//                if (haveRegion.insert(region))
//                {
//                    nOriginating[proci]++;
//                }
//            }
//        }
//    }
//
//    if (debug)
//    {
//        Pout<< "Counted " << nOriginating[Pstream::myProcNo()]
//            << " local regions." << endl;
//    }
//
//
//    // Global numbering for compacted local regions
//    autoPtr<globalIndex> globalCompactPtr
//    (
//        new globalIndex(nOriginating[Pstream::myProcNo()])
//    );
//    const globalIndex& globalCompact = globalCompactPtr();
//
//
//    // 4b: renumber
//    // Renumber into compact indices. Note that since we've already made
//    // all regions global we now need a Map to store the compacting
//    // information
//    // instead of a labelList - otherwise we could have used a straight
//    // labelList.
//
//    // Local compaction map
//    Map<label> globalToCompact(2*nOriginating[Pstream::myProcNo()]);
//    // Remote regions we want the compact number for
//    List<labelHashSet> nonLocal(Pstream::nProcs());
//    forAll(nonLocal, proci)
//    {
//        if (proci != Pstream::myProcNo())
//        {
//            nonLocal[proci].resize(2*nOriginating[proci]);
//        }
//    }
//
//    forAll(cellRegion, celli)
//    {
//        label region = cellRegion[celli];
//        if (globalRegions.isLocal(region))
//        {
//            // Insert new compact region (if not yet present)
//            globalToCompact.insert
//            (
//                region,
//                globalCompact.toGlobal(globalToCompact.size())
//            );
//        }
//        else
//        {
//            nonLocal[globalRegions.whichProcID(region)].insert(region);
//        }
//    }
//
//
//    // Now we have all the local regions compacted. Now we need to get the
//    // non-local ones from the processors to whom they are local.
//    // Convert the nonLocal (labelHashSets) to labelLists.
//
//    labelListList sendNonLocal(Pstream::nProcs());
//    forAll(sendNonLocal, proci)
//    {
//        sendNonLocal[proci] = nonLocal[proci].toc();
//    }
//
//    if (debug)
//    {
//        forAll(sendNonLocal, proci)
//        {
//            Pout<< "    from processor " << proci
//                << " want " << sendNonLocal[proci].size()
//                << " region numbers." << endl;
//        }
//        Pout<< endl;
//    }
//
//
//    // Get the wanted region labels into recvNonLocal
//    labelListList recvNonLocal;
//    Pstream::exchange<labelList, label>(sendNonLocal, recvNonLocal);
//
//    // Now we have the wanted compact region labels that proci wants in
//    // recvNonLocal[proci]. Construct corresponding list of compact
//    // region labels to send back.
//
//    labelListList sendWantedLocal(Pstream::nProcs());
//    forAll(recvNonLocal, proci)
//    {
//        const labelList& nonLocal = recvNonLocal[proci];
//        sendWantedLocal[proci].setSize(nonLocal.size());
//
//        forAll(nonLocal, i)
//        {
//            sendWantedLocal[proci][i] = globalToCompact[nonLocal[i]];
//        }
//    }
//
//
//    // Send back (into recvNonLocal)
//    recvNonLocal.clear();
//    Pstream::exchange<labelList, label>(sendWantedLocal, recvNonLocal);
//    sendWantedLocal.clear();
//
//    // Now recvNonLocal contains for every element in setNonLocal the
//    // corresponding compact number. Insert these into the local compaction
//    // map.
//
//    forAll(recvNonLocal, proci)
//    {
//        const labelList& wantedRegions = sendNonLocal[proci];
//        const labelList& compactRegions = recvNonLocal[proci];
//
//        forAll(wantedRegions, i)
//        {
//            globalToCompact.insert(wantedRegions[i], compactRegions[i]);
//        }
//    }
//
//    // Finally renumber the regions
//    forAll(cellRegion, celli)
//    {
//        cellRegion[celli] = globalToCompact[cellRegion[celli]];
//    }
//
//    return globalCompactPtr;
//}


void Foam::cellCellStencils::inverseDistance::findHoles
(
    const globalIndex& globalCells,
    const fvMesh& mesh,
    const labelList& zoneID,
    const labelListList& stencil,
    labelList& cellTypes
) const
{
    const fvBoundaryMesh& pbm = mesh.boundary();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();


    // The input cellTypes will be
    // - HOLE           : cell part covered by other-mesh patch
    // - INTERPOLATED   : cell fully covered by other-mesh patch
    //                    or next to 'overset' patch
    // - CALCULATED     : otherwise
    //
    // so we start a walk from our patches and any cell we cannot reach
    // (because we walk is stopped by other-mesh patch) is a hole.


    DebugInfo<< FUNCTION_NAME << " : Starting hole flood filling" << endl;

    DebugInfo<< FUNCTION_NAME << " : Starting hole cells : "
        << findIndices(cellTypes, HOLE).size() << endl;

    boolList isBlockedFace(mesh.nFaces(), false);
    label nBlocked = 0;

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label ownType = cellTypes[own[faceI]];
        label neiType = cellTypes[nei[faceI]];
        if
        (
             (ownType == HOLE && neiType != HOLE)
          || (ownType != HOLE && neiType == HOLE)
        )
        {
            isBlockedFace[faceI] = true;
            nBlocked++;
        }
    }
    DebugInfo<< FUNCTION_NAME << " : Marked internal hole boundaries : "
        << nBlocked << endl;


    labelList nbrCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, nbrCellTypes);

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label ownType = cellTypes[own[faceI]];
        label neiType = nbrCellTypes[faceI-mesh.nInternalFaces()];

        if
        (
             (ownType == HOLE && neiType != HOLE)
          || (ownType != HOLE && neiType == HOLE)
        )
        {
            isBlockedFace[faceI] = true;
            nBlocked++;
        }
    }

    DebugInfo<< FUNCTION_NAME << " : Marked all hole boundaries : "
        << nBlocked << endl;

    // Determine regions
    regionSplit cellRegion(mesh, isBlockedFace);
    const label nRegions = cellRegion.nRegions();

    //labelList cellRegion;
    //label nRegions = -1;
    //{
    //    const globalIndex globalFaces(mesh.nFaces());
    //    uncompactedRegionSplit
    //    (
    //        mesh,
    //        globalFaces,
    //        gMax(zoneID)+1,
    //        zoneID,
    //        cellTypes,
    //        isBlockedFace,
    //        cellRegion
    //    );
    //    autoPtr<globalIndex> globalRegions
    //    (
    //        compactedRegionSplit
    //        (
    //            mesh,
    //            globalFaces,
    //            cellRegion
    //        )
    //    );
    //    nRegions = globalRegions().size();
    //}
    DebugInfo<< FUNCTION_NAME << " : Determined regions : "
        << nRegions << endl;

    //Info<< typeName << " : detected " << nRegions
    //    << " mesh regions after overset" << nl << endl;



    // Now we'll have a mesh split according to where there are cells
    // covered by the other-side patches. See what we can reach from our
    // real patches

    //  0 : region not yet determined
    //  1 : borders blockage so is not ok (but can be overridden by real
    //      patch)
    //  2 : has real patch in it so is reachable
    labelList regionType(nRegions, Zero);


    // See if any regions borders blockage. Note: isBlockedFace is already
    // parallel synchronised.
    {
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            if (isBlockedFace[faceI])
            {
                label ownRegion = cellRegion[own[faceI]];

                if (cellTypes[own[faceI]] != HOLE)
                {
                    if (regionType[ownRegion] == 0)
                    {
                        regionType[ownRegion] = 1;
                    }
                }

                label neiRegion = cellRegion[nei[faceI]];

                if (cellTypes[nei[faceI]] != HOLE)
                {
                    if (regionType[neiRegion] == 0)
                    {
                        regionType[neiRegion] = 1;
                    }
                }
            }
        }
        for
        (
            label faceI = mesh.nInternalFaces();
            faceI < mesh.nFaces();
            faceI++
        )
        {
            if (isBlockedFace[faceI])
            {
                label ownRegion = cellRegion[own[faceI]];

                if (regionType[ownRegion] == 0)
                {
                    regionType[ownRegion] = 1;
                }
            }
        }
    }


    // Override with real patches
    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];

        if (isA<oversetFvPatch>(fvp))
        {}
        else if (!fvPatch::constraintType(fvp.type()))
        {
            const labelList& fc = fvp.faceCells();
            forAll(fc, i)
            {
                label regionI = cellRegion[fc[i]];

                if (cellTypes[fc[i]] != HOLE && regionType[regionI] != 2)
                {
                    regionType[regionI] = 2;
                }
            }
        }
    }

    DebugInfo<< FUNCTION_NAME << " : Done local analysis" << endl;

    // Now we've handled
    // - cells next to blocked cells
    // - coupled boundaries
    // Only thing to handle is the interpolation between regions


    labelListList compactStencil(stencil);
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, compactStencil, compactMap);

    DebugInfo<< FUNCTION_NAME << " : Converted stencil into compact form"
        << endl;


    while (true)
    {
        // Synchronise region status on processors
        // (could instead swap status through processor patches)
        Pstream::listCombineGather(regionType, maxEqOp<label>());
        Pstream::listCombineScatter(regionType);

        DebugInfo<< FUNCTION_NAME << " : Gathered region type" << endl;

        // Communicate region status through interpolative cells
        labelList cellRegionType(UIndirectList<label>(regionType, cellRegion));
        map.distribute(cellRegionType);

        DebugInfo<< FUNCTION_NAME << " : Interpolated region type" << endl;



        label nChanged = 0;
        forAll(pbm, patchI)
        {
            const fvPatch& fvp = pbm[patchI];

            if (isA<oversetFvPatch>(fvp))
            {
                const labelUList& fc = fvp.faceCells();
                forAll(fc, i)
                {
                    label cellI = fc[i];
                    label regionI = cellRegion[cellI];

                    if (regionType[regionI] != 2)
                    {
                        const labelList& slots = compactStencil[cellI];
                        forAll(slots, i)
                        {
                            label otherType = cellRegionType[slots[i]];

                            if (otherType == 2)
                            {
                                //Pout<< "Reachable through interpolation : "
                                //    << regionI << " at cell "
                                //    << mesh.cellCentres()[cellI] << endl;
                                regionType[regionI] = 2;
                                nChanged++;
                                break;
                            }
                        }
                    }
                }
            }
        }

        reduce(nChanged, sumOp<label>());
        DebugInfo<< FUNCTION_NAME << " : Determined regions changed : "
            << nChanged << endl;

        if (nChanged == 0)
        {
            break;
        }
    }


    // See which regions have not been visited (regionType == 1)
    forAll(cellRegion, cellI)
    {
        label type = regionType[cellRegion[cellI]];
        if (type == 1 && cellTypes[cellI] != HOLE)
        {
            cellTypes[cellI] = HOLE;
        }
    }
    DebugInfo<< FUNCTION_NAME << " : Finished hole flood filling" << endl;
}


void Foam::cellCellStencils::inverseDistance::seedCell
(
    const label cellI,
    const scalar wantedFraction,
    PackedBoolList& isFront,
    scalarField& fraction
) const
{
    const cell& cFaces = mesh_.cells()[cellI];
    forAll(cFaces, i)
    {
        label nbrFacei = cFaces[i];
        if (fraction[nbrFacei] < wantedFraction)
        {
            fraction[nbrFacei] = wantedFraction;
            isFront.set(nbrFacei);
        }
    }
}


void Foam::cellCellStencils::inverseDistance::walkFront
(
    const scalar layerRelax,
    const labelListList& allStencil,
    labelList& allCellTypes,
    scalarField& allWeight
) const
{
    // Current front
    PackedBoolList isFront(mesh_.nFaces());

    const fvBoundaryMesh& fvm = mesh_.boundary();


    // 'overset' patches

    forAll(fvm, patchI)
    {
        if (isA<oversetFvPatch>(fvm[patchI]))
        {
            const labelList& fc = fvm[patchI].faceCells();
            forAll(fc, i)
            {
                label cellI = fc[i];
                if (allCellTypes[cellI] == INTERPOLATED)
                {
                    // Note that acceptors might have been marked hole if
                    // there are no donors in which case we do not want to
                    // walk this out. This is an extreme situation.
                    isFront.set(fvm[patchI].start()+i);
                }
            }
        }
    }


    // Outside of 'hole' region
    {
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label ownType = allCellTypes[own[faceI]];
            label neiType = allCellTypes[nei[faceI]];
            if
            (
                 (ownType == HOLE && neiType != HOLE)
              || (ownType != HOLE && neiType == HOLE)
            )
            {
                //Pout<< "Front at face:" << faceI
                //    << " at:" << mesh_.faceCentres()[faceI] << endl;
                isFront.set(faceI);
            }
        }

        labelList nbrCellTypes;
        syncTools::swapBoundaryCellList(mesh_, allCellTypes, nbrCellTypes);

        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            label ownType = allCellTypes[own[faceI]];
            label neiType = nbrCellTypes[faceI-mesh_.nInternalFaces()];

            if
            (
                 (ownType == HOLE && neiType != HOLE)
              || (ownType != HOLE && neiType == HOLE)
            )
            {
                //Pout<< "Front at coupled face:" << faceI
                //    << " at:" << mesh_.faceCentres()[faceI] << endl;
                isFront.set(faceI);
            }
        }
    }


    // Current interpolation fraction
    scalarField fraction(mesh_.nFaces(), Zero);

    forAll(isFront, faceI)
    {
        if (isFront.get(faceI))
        {
            fraction[faceI] = 1.0;
        }
    }


    while (returnReduce(isFront.count() > 0, orOp<bool>()))
    {
        // Interpolate cells on front
        PackedBoolList newIsFront(mesh_.nFaces());
        scalarField newFraction(fraction);
        forAll(isFront, faceI)
        {
            if (isFront.get(faceI))
            {
                label own = mesh_.faceOwner()[faceI];
                if (allCellTypes[own] != HOLE)
                {
                    if (allWeight[own] < fraction[faceI])
                    {
                        // Cell wants to become interpolated (if sufficient
                        // stencil, otherwise becomes hole)
                        if (allStencil[own].size())
                        {
                            allWeight[own] = fraction[faceI];
                            allCellTypes[own] = INTERPOLATED;
                            // Add faces of cell (with lower weight) as new
                            // front
                            seedCell
                            (
                                own,
                                fraction[faceI]-layerRelax,
                                newIsFront,
                                newFraction
                            );
                        }
                        else
                        {
                            allWeight[own] = 0.0;
                            allCellTypes[own] = HOLE;
                            // Add faces of cell as new front
                            seedCell
                            (
                                own,
                                1.0,
                                newIsFront,
                                newFraction
                            );
                        }
                    }
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (allCellTypes[nei] != HOLE)
                    {
                        if (allWeight[nei] < fraction[faceI])
                        {
                            if (allStencil[nei].size())
                            {
                                allWeight[nei] = fraction[faceI];
                                allCellTypes[nei] = INTERPOLATED;
                                seedCell
                                (
                                    nei,
                                    fraction[faceI]-layerRelax,
                                    newIsFront,
                                    newFraction
                                );
                            }
                            else
                            {
                                allWeight[nei] = 0.0;
                                allCellTypes[nei] = HOLE;
                                seedCell
                                (
                                    nei,
                                    1.0,
                                    newIsFront,
                                    newFraction
                                );
                            }
                        }
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh_, newIsFront, orEqOp<unsigned int>());
        syncTools::syncFaceList(mesh_, newFraction, maxEqOp<scalar>());

        isFront.transfer(newIsFront);
        fraction.transfer(newFraction);
    }
}


void Foam::cellCellStencils::inverseDistance::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    // Inverse-distance weighting

    weights.setSize(donorCcs.size());
    scalar sum = 0.0;
    forAll(donorCcs, i)
    {
        scalar d = mag(sample-donorCcs[i]);

        if (d > ROOTVSMALL)
        {
            weights[i] = 1.0/d;
            sum += weights[i];
        }
        else
        {
            // Short circuit
            weights = 0.0;
            weights[i] = 1.0;
            return;
        }
    }
    forAll(weights, i)
    {
        weights[i] /= sum;
    }
}


void Foam::cellCellStencils::inverseDistance::createStencil
(
    const globalIndex& globalCells
)
{
    // Send cell centre back to donor
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // The complication is that multiple acceptors need the same donor
    // (but with different weights obviously)
    // So we do multi-pass:
    // - send over cc of acceptor for which we want stencil.
    //   Consistently choose the acceptor with smallest magSqr in case of
    //   multiple acceptors for the containing cell/donor.
    // - find the cell-cells and weights for the donor
    // - send back together with the acceptor cc
    // - use the acceptor cc to see if it was 'me' that sent it. If so
    //   mark me as complete so it doesn't get involved in the next loop.
    // - loop until all resolved.

    // Special value for unused points
    const vector greatPoint(GREAT, GREAT, GREAT);

    boolList isValidDonor(mesh_.nCells(), true);
    forAll(cellTypes_, celli)
    {
        if (cellTypes_[celli] == HOLE)
        {
            isValidDonor[celli] = false;
        }
    }


    // Has acceptor been handled already?
    PackedBoolList doneAcceptor(interpolationCells_.size());

    while (true)
    {
        pointField samples(cellInterpolationMap().constructSize(), greatPoint);

        // Fill remote slots (override old content). We'll find out later
        // on which one has won and mark this one in doneAcceptor.
        label nSamples = 0;
        forAll(interpolationCells_, i)
        {
            if (!doneAcceptor[i])
            {
                label cellI = interpolationCells_[i];
                const point& cc = mesh_.cellCentres()[cellI];
                const labelList& slots = cellStencil_[cellI];

                if (slots.size() != 1)
                {
                    FatalErrorInFunction<< "Problem:" << slots
                        << abort(FatalError);
                }

                forAll(slots, slotI)
                {
                    label elemI = slots[slotI];
                    //Pout<< "    acceptor:" << cellI
                    //    << " at:" << mesh_.cellCentres()[cellI]
                    //    << " global:" << globalCells.toGlobal(cellI)
                    //    << " found in donor:" << elemI << endl;
                    minMagSqrEqOp<point>()(samples[elemI], cc);
                }
                nSamples++;
            }
        }


        if (returnReduce(nSamples, sumOp<label>()) == 0)
        {
            break;
        }

        // Send back to donor. Make sure valid point takes priority
        mapDistributeBase::distribute<point, minMagSqrEqOp<point>, flipOp>
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            mesh_.nCells(),
            cellInterpolationMap().constructMap(),
            false,
            cellInterpolationMap().subMap(),
            false,
            samples,
            minMagSqrEqOp<point>(),
            flipOp(),                               // negateOp
            greatPoint,                             // nullValue
            UPstream::msgType()//,
//             cellInterpolationMap().comm()
        );

        // All the donor cells will now have a valid cell centre. Construct a
        // stencil for these.

        DynamicList<label> donorCells(mesh_.nCells());
        forAll(samples, cellI)
        {
            if (samples[cellI] != greatPoint)
            {
                donorCells.append(cellI);
            }
        }


        // Get neighbours (global cell and centre) of donorCells.
        labelListList donorCellCells(mesh_.nCells());
        pointListList donorCellCentres(mesh_.nCells());
        globalCellCells
        (
            globalCells,
            mesh_,
            isValidDonor,
            donorCells,
            donorCellCells,
            donorCellCentres
        );

        // Determine the weights.
        scalarListList donorWeights(mesh_.nCells());
        forAll(donorCells, i)
        {
            label cellI = donorCells[i];
            const pointList& donorCentres = donorCellCentres[cellI];
            stencilWeights
            (
                samples[cellI],
                donorCentres,
                donorWeights[cellI]
            );
        }

        // Transfer the information back to the acceptor:
        // - donorCellCells : stencil (with first element the original donor)
        // - donorWeights : weights for donorCellCells
        cellInterpolationMap().distribute(donorCellCells);
        cellInterpolationMap().distribute(donorWeights);
        cellInterpolationMap().distribute(samples);

        // Check which acceptor has won and transfer
        forAll(interpolationCells_, i)
        {
            if (!doneAcceptor[i])
            {
                label cellI = interpolationCells_[i];
                const labelList& slots = cellStencil_[cellI];

                if (slots.size() != 1)
                {
                    FatalErrorInFunction << "Problem:" << slots
                        << abort(FatalError);
                }

                label slotI = slots[0];

                // Important: check if the stencil is actually for this cell
                if (samples[slotI] == mesh_.cellCentres()[cellI])
                {
                    cellStencil_[cellI].transfer(donorCellCells[slotI]);
                    cellInterpolationWeights_[cellI].transfer
                    (
                        donorWeights[slotI]
                    );
                    // Mark cell as being done so it does not get sent over
                    // again.
                    doneAcceptor.set(i);
                }
            }
        }
    }

    // Re-do the mapDistribute
    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute
        (
            globalCells,
            cellStencil_,
            compactMap
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::inverseDistance::inverseDistance
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    cellCellStencil(mesh),
    dict_(dict),
    smallVec_(Zero),
    cellTypes_(labelList(mesh.nCells(), CALCULATED)),
    interpolationCells_(0),
    cellInterpolationMap_(),
    cellStencil_(0),
    cellInterpolationWeights_(0),
    cellInterpolationWeight_
    (
        IOobject
        (
            "cellInterpolationWeight",
            mesh_.facesInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{
    // Protect local fields from interpolation
    nonInterpolatedFields_.insert("cellInterpolationWeight");
    nonInterpolatedFields_.insert("cellTypes");
    nonInterpolatedFields_.insert("maxMagWeight");

    // For convenience also suppress frequently used displacement field
    nonInterpolatedFields_.insert("cellDisplacement");
    nonInterpolatedFields_.insert("grad(cellDisplacement)");
    const word w("snGradCorr(cellDisplacement)");
    const word d("((viscosity*faceDiffusivity)*magSf)");
    nonInterpolatedFields_.insert("surfaceIntegrate(("+d+"*"+w+"))");

    // Read zoneID
    this->zoneID();

    // Read old-time cellTypes
    IOobject io
    (
        "cellTypes",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );
    if (io.typeHeaderOk<volScalarField>(true))
    {
        if (debug)
        {
            Pout<< "Reading cellTypes from time " << mesh_.time().timeName()
                << endl;
        }

        const volScalarField volCellTypes(io, mesh_);
        forAll(volCellTypes, celli)
        {
            // Round to integer
            cellTypes_[celli] = volCellTypes[celli];
        }
    }

    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::inverseDistance::~inverseDistance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCellStencils::inverseDistance::update()
{
    scalar layerRelax(dict_.lookupOrDefault("layerRelax", 1.0));

    scalar tol = dict_.lookupOrDefault("tolerance", 1e-10);
    smallVec_ = mesh_.bounds().span()*tol;

    const labelIOList& zoneID = this->zoneID();

    label nZones = gMax(zoneID)+1;
    labelList nCellsPerZone(nZones, Zero);
    forAll(zoneID, cellI)
    {
        nCellsPerZone[zoneID[cellI]]++;
    }
    Pstream::listCombineGather(nCellsPerZone, plusEqOp<label>());
    Pstream::listCombineScatter(nCellsPerZone);


    const boundBox& allBb(mesh_.bounds());


    PtrList<fvMeshSubset> meshParts(nZones);
    List<treeBoundBoxList> meshBb(nZones);

    // Determine zone meshes and bounding boxes
    {
        // Per processor, per zone the bounding box
        List<treeBoundBoxList> procBb(Pstream::nProcs());
        procBb[Pstream::myProcNo()].setSize(nZones);

        forAll(meshParts, zonei)
        {
            meshParts.set
            (
                zonei,
                new fvMeshSubset(mesh_)
            );
            meshParts[zonei].setLargeCellSubset(zoneID, zonei);
            const fvMesh& subMesh = meshParts[zonei].subMesh();

            // Trigger early evaluation of mesh dimension (in case there are
            // zero cells in mesh)
            (void)subMesh.nGeometricD();

            if (subMesh.nPoints())
            {
                procBb[Pstream::myProcNo()][zonei] =
                    treeBoundBox(subMesh.points());
                procBb[Pstream::myProcNo()][zonei].inflate(1e-6);
            }
            else
            {
                // No part of zone on this processor. Make up bb.
                procBb[Pstream::myProcNo()][zonei] = treeBoundBox
                (
                    allBb.min() - 2*allBb.span(),
                    allBb.min() - allBb.span()
                );
                procBb[Pstream::myProcNo()][zonei].inflate(1e-6);
            }
        }

        Pstream::gatherList(procBb);
        Pstream::scatterList(procBb);

        // Move local bounding boxes to per-mesh indexing
        forAll(meshBb, zoneI)
        {
            treeBoundBoxList& bbs = meshBb[zoneI];
            bbs.setSize(Pstream::nProcs());
            forAll(procBb, procI)
            {
                bbs[procI] = procBb[procI][zoneI];
            }
        }
    }


    // Determine patch bounding boxes. These are either global and provided
    // by the user or processor-local as a copy of the mesh bounding box.

    List<treeBoundBoxList> patchBb(nZones);
    List<labelVector> patchDivisions(nZones);
    PtrList<PackedList<2>> patchParts(nZones);
    labelList allPatchTypes(mesh_.nCells(), OTHER);

    {
        treeBoundBox globalPatchBb;
        if (dict_.readIfPresent("searchBox", globalPatchBb))
        {
            // All processors, all zones have the same bounding box
            patchBb = treeBoundBoxList(Pstream::nProcs(), globalPatchBb);
        }
        else
        {
            // Use the meshBb (differing per zone, per processor)
            patchBb = meshBb;
        }
    }

    {
        labelVector globalDivs;
        if (dict_.readIfPresent("searchBoxDivisions", globalDivs))
        {
            patchDivisions = globalDivs;
        }
        else
        {
            const labelVector& dim = mesh_.geometricD();
            label nDivs = -1;
            if (mesh_.nGeometricD() == 1)
            {
                nDivs = mesh_.nCells();
            }
            else if (mesh_.nGeometricD() == 2)
            {
                nDivs = label(Foam::sqrt(scalar(mesh_.nCells())));
            }
            else
            {
                nDivs = label(Foam::cbrt(scalar(mesh_.nCells())));
            }

            labelVector v(nDivs, nDivs, nDivs);
            forAll(dim, i)
            {
                if (dim[i] == -1)
                {
                    v[i] = 1;
                }
            }
            patchDivisions = v;
        }
    }

    forAll(patchParts, zoneI)
    {
        patchParts.set
        (
            zoneI,
            new PackedList<2>
            (
                patchDivisions[zoneI][0]
               *patchDivisions[zoneI][1]
               *patchDivisions[zoneI][2]
            )
        );
        markBoundaries
        (
            meshParts[zoneI].subMesh(),
            smallVec_,

            patchBb[zoneI][Pstream::myProcNo()],
            patchDivisions[zoneI],
            patchParts[zoneI],

            meshParts[zoneI].cellMap(),
            allPatchTypes
        );
    }


    // Print a bit
    if (debug)
    {
        Info<< type() << " : detected " << nZones
            << " mesh regions" << endl;
        Info<< incrIndent;
        forAll(nCellsPerZone, zoneI)
        {
            Info<< indent<< "zone:" << zoneI
                << " nCells:" << nCellsPerZone[zoneI]
                << "  voxels:" << patchDivisions[zoneI]
                << " bb:" << patchBb[zoneI][Pstream::myProcNo()]
                << endl;
        }
        Info<< decrIndent;
    }


    // Current best guess for cells. Includes best stencil. Weights should
    // add up to volume.
    labelList allCellTypes(mesh_.nCells(), CALCULATED);
    labelListList allStencil(mesh_.nCells());
    // zoneID of donor
    labelList allDonorID(mesh_.nCells(), -1);

    const globalIndex globalCells(mesh_.nCells());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Mark holes (in allCellTypes)
    for (label srcI = 0; srcI < meshParts.size()-1; srcI++)
    {
        for (label tgtI = srcI+1; tgtI < meshParts.size(); tgtI++)
        {
            markPatchesAsHoles
            (
                pBufs,

                meshParts,

                patchBb,
                patchDivisions,
                patchParts,

                srcI,
                tgtI,
                allCellTypes
            );
            markPatchesAsHoles
            (
                pBufs,

                meshParts,

                patchBb,
                patchDivisions,
                patchParts,

                tgtI,
                srcI,
                allCellTypes
            );
        }
    }

    // Find donors (which are not holes) in allStencil, allDonorID
    for (label srcI = 0; srcI < meshParts.size()-1; srcI++)
    {
        for (label tgtI = srcI+1; tgtI < meshParts.size(); tgtI++)
        {
            markDonors
            (
                globalCells,
                pBufs,
                meshParts,
                meshBb,
                allCellTypes,

                tgtI,
                srcI,
                allStencil,
                allDonorID
            );
            markDonors
            (
                globalCells,
                pBufs,
                meshParts,
                meshBb,
                allCellTypes,

                srcI,
                tgtI,
                allStencil,
                allDonorID
            );
        }
    }

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes", allCellTypes)
        );
        tfld().write();
    }

    // Use the patch types and weights to decide what to do
    forAll(allPatchTypes, cellI)
    {
        if (allCellTypes[cellI] != HOLE)
        {
            switch (allPatchTypes[cellI])
            {
                case OVERSET:
                {
                    // Require interpolation. See if possible.

                    if (allStencil[cellI].size())
                    {
                        allCellTypes[cellI] = INTERPOLATED;
                    }
                    else
                    {
                        // If there are no donors we can either set the
                        // acceptor cell to 'hole' i.e. nothing gets calculated
                        // if the acceptor cells go outside the donor meshes or
                        // we can set it to 'calculated' to have something
                        // like an acmi type behaviour where only covered
                        // acceptor cell use the interpolation and non-covered
                        // become calculated. Uncomment below line. Note: this
                        // should be switchable?
                        //allCellTypes[cellI] = CALCULATED;

                        allCellTypes[cellI] = HOLE;
                    }
                }
            }
        }
    }

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_patch", allCellTypes)
        );
        tfld().write();
    }


    // Check previous iteration cellTypes_ for any hole->calculated changes
    // If so set the cell either to interpolated (if there are donors) or
    // holes (if there are no donors). Note that any interpolated cell might
    // still be overwritten by the flood filling
    {
        label nCalculated = 0;

        forAll(cellTypes_, celli)
        {
            if (allCellTypes[celli] == CALCULATED && cellTypes_[celli] == HOLE)
            {
                if (allStencil[celli].size() == 0)
                {
                    // Reset to hole
                    allCellTypes[celli] = HOLE;
                    allStencil[celli].clear();
                }
                else
                {
                    allCellTypes[celli] = INTERPOLATED;
                    nCalculated++;
                }
            }
        }

        if (debug)
        {
            Pout<< "Detected " << nCalculated << " cells changing from hole"
                << " to calculated. Changed to interpolated"
                << endl;
        }
    }


    // Mark unreachable bits
    findHoles(globalCells, mesh_, zoneID, allStencil, allCellTypes);

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_hole", allCellTypes)
        );
        tfld().write();
    }
    if (debug)
    {
        labelList stencilSize(mesh_.nCells());
        forAll(allStencil, celli)
        {
            stencilSize[celli] = allStencil[celli].size();
        }
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allStencil_hole", stencilSize)
        );
        tfld().write();
    }


    // Add buffer interpolation layer(s) around holes
    scalarField allWeight(mesh_.nCells(), Zero);
    walkFront(layerRelax, allStencil, allCellTypes, allWeight);

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_front", allCellTypes)
        );
        tfld().write();
    }


    // Convert cell-cell addressing to stencil in compact notation

    cellTypes_.transfer(allCellTypes);
    cellStencil_.setSize(mesh_.nCells());
    cellInterpolationWeights_.setSize(mesh_.nCells());
    DynamicList<label> interpolationCells;
    forAll(cellTypes_, cellI)
    {
        if (cellTypes_[cellI] == INTERPOLATED)
        {
            cellStencil_[cellI].transfer(allStencil[cellI]);
            cellInterpolationWeights_[cellI].setSize(1);
            cellInterpolationWeights_[cellI][0] = 1.0;
            interpolationCells.append(cellI);
        }
        else
        {
            cellStencil_[cellI].clear();
            cellInterpolationWeights_[cellI].clear();
        }
    }
    interpolationCells_.transfer(interpolationCells);

    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute(globalCells, cellStencil_, compactMap)
    );
    cellInterpolationWeight_.transfer(allWeight);
    dynamicOversetBlastFvMesh::correctBoundaryConditions
    <
        volScalarField,
        oversetFvPatchField<scalar>
    >(cellInterpolationWeight_.boundaryFieldRef(), false);


    if (debug&2)
    {
        // Dump mesh
        mesh_.time().write();

        // Dump stencil
        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"injectionStencil.obj");
        Pout<< type() << " : dumping injectionStencil to "
            << str.name() << endl;
        pointField cc(mesh_.cellCentres());
        cellInterpolationMap().distribute(cc);

        forAll(cellStencil_, celli)
        {
            const labelList& slots = cellStencil_[celli];
            if (slots.size())
            {
                const point& accCc = mesh_.cellCentres()[celli];
                forAll(slots, i)
                {
                    const point& donorCc = cc[slots[i]];
                    str.write(linePointRef(accCc, 0.1*accCc+0.9*donorCc));
                }
            }
        }
    }


    // Extend stencil to get inverse distance weighted neighbours
    createStencil(globalCells);


    if (debug&2)
    {
        // Dump weight
        cellInterpolationWeight_.instance() = mesh_.time().timeName();
        cellInterpolationWeight_.write();

        // Dump max weight
        {
            scalarField maxMagWeight(mesh_.nCells(), Zero);
            forAll(cellStencil_, celli)
            {
                const scalarList& wghts = cellInterpolationWeights_[celli];
                forAll(wghts, i)
                {
                    if (mag(wghts[i]) > mag(maxMagWeight[celli]))
                    {
                        maxMagWeight[celli] = wghts[i];
                    }
                }
                if (mag(maxMagWeight[celli]) > 1)
                {
                    const pointField& cc = mesh_.cellCentres();
                    Pout<< "cell:" << celli
                        << " at:" << cc[celli]
                        << " zone:" << zoneID[celli]
                        << " donors:" << cellStencil_[celli]
                        << " weights:" << wghts
                        << " coords:"
                        << UIndirectList<point>(cc, cellStencil_[celli])
                        << " donorZone:"
                        << UIndirectList<label>(zoneID, cellStencil_[celli])
                        << endl;
                }
            }
            tmp<volScalarField> tfld
            (
                createField(mesh_, "maxMagWeight", maxMagWeight)
            );
            dynamicOversetBlastFvMesh::correctBoundaryConditions
            <
                volScalarField,
                oversetFvPatchField<scalar>
            >(tfld.ref().boundaryFieldRef(), false);
            tfld().write();
        }

        // Dump cell types
        {
            tmp<volScalarField> tfld
            (
                createField(mesh_, "cellTypes", cellTypes_)
            );
            //tfld.ref().correctBoundaryConditions();
            dynamicOversetBlastFvMesh::correctBoundaryConditions
            <
                volScalarField,
                oversetFvPatchField<scalar>
            >(tfld.ref().boundaryFieldRef(), false);
            tfld().write();
        }


        // Dump stencil, one per zone
        mkDir(mesh_.time().timePath());
        pointField cc(mesh_.cellCentres());
        cellInterpolationMap().distribute(cc);
        forAll(meshParts, zonei)
        {
            OBJstream str
            (
               mesh_.time().timePath()
               + "/stencil_" + name(zonei) + ".obj"
            );
            Pout<< type() << " : dumping to " << str.name() << endl;

            const labelList& subMeshCellMap = meshParts[zonei].cellMap();

            forAll(subMeshCellMap, subcelli)
            {
                const label celli = subMeshCellMap[subcelli];
                const labelList& slots = cellStencil_[celli];
                const point& accCc = mesh_.cellCentres()[celli];
                forAll(slots, i)
                {
                    const point& donorCc = cc[slots[i]];
                    str.write(linePointRef(accCc, 0.1*accCc+0.9*donorCc));
                }
            }
        }
    }

    // Print some stats
    {
        labelList nCells(count(3, cellTypes_));

        label nLocal = 0;
        label nMixed = 0;
        label nRemote = 0;
        forAll(interpolationCells_, i)
        {
            label celli = interpolationCells_[i];
            const labelList& slots = cellStencil_[celli];

            bool hasLocal = false;
            bool hasRemote = false;

            forAll(slots, sloti)
            {
                if (slots[sloti] >= mesh_.nCells())
                {
                    hasRemote = true;
                }
                else
                {
                    hasLocal = true;
                }
            }

            if (hasRemote)
            {
                if (!hasLocal)
                {
                    nRemote++;
                }
                else
                {
                    nMixed++;
                }
            }
            else if (hasLocal)
            {
                nLocal++;
            }
        }
        reduce(nLocal, sumOp<label>());
        reduce(nMixed, sumOp<label>());
        reduce(nRemote, sumOp<label>());

        Info<< "Overset analysis : nCells : "
            << returnReduce(cellTypes_.size(), sumOp<label>()) << nl
            << incrIndent
            << indent << "calculated   : " << nCells[CALCULATED] << nl
            << indent << "interpolated : " << nCells[INTERPOLATED]
            << " (interpolated from local:" << nLocal
            << "  mixed local/remote:" << nMixed
            << "  remote:" << nRemote << ")" << nl
            << indent << "hole         : " << nCells[HOLE] << nl
            << decrIndent << endl;
    }

    // Tbd: detect if anything changed. Most likely it did!
    return true;
}


// ************************************************************************* //
