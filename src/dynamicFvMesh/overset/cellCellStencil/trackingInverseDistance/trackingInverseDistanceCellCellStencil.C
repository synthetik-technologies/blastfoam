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

#include "trackingInverseDistanceCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
#include "Time.H"
#include "fvMeshSubset.H"
#include "globalIndex.H"
#include "oversetFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "treeBoundBoxList.H"
#include "voxelMeshSearch.H"
#include "dynamicOversetBlastFvMesh.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(trackingInverseDistance, 0);
    addToRunTimeSelectionTable(cellCellStencil, trackingInverseDistance, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cellCellStencils::trackingInverseDistance::markBoundaries
(
    const fvMesh& mesh,
    const vector& smallVec,

    const boundBox& bb,
    labelVector& nDivs,
    PackedList<2>& patchTypes,

    const labelList& cellMap,
    labelList& patchCellTypes
)
{
    // Mark all voxels that overlap the bounding box of any patch

    static bool hasWarned = false;

    const fvBoundaryMesh& pbm = mesh.boundary();

    patchTypes = patchCellType::OTHER;

    // Mark wall boundaries
    forAll(pbm, patchi)
    {
        const fvPatch& fvp = pbm[patchi];
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
                    voxelMeshSearch::fill
                    (
                        patchTypes,
                        bb,
                        nDivs,
                        faceBb,
                        patchCellType::PATCH
                    );
                }
            }
        }
    }

    // Override with overset boundaries
    forAll(pbm, patchi)
    {
        const fvPatch& fvp = pbm[patchi];
        const labelList& fc = fvp.faceCells();

        if (isA<oversetFvPatch>(fvp))
        {
            //Info<< "Marking cells on overset patch " << fvp.name() << endl;
            const polyPatch& pp = fvp.patch();
            const vectorField::subField faceCentres(pp.faceCentres());
            forAll(pp, i)
            {
                // Mark in overall patch types
                patchCellTypes[cellMap[fc[i]]] = patchCellType::OVERSET;

                // Mark in voxel mesh
                boundBox faceBb(pp.points(), pp[i], false);
                faceBb.min() -= smallVec;
                faceBb.max() += smallVec;
                if (!bb.contains(faceCentres[i]))
                {
                    if (!hasWarned)
                    {
                        hasWarned = true;
                        WarningInFunction << "Patch " << pp.name()
                            << " : face at " << faceCentres[i]
                            << " is not inside search bounding box of"
                            << " voxel mesh " << bb << endl
                            << "    Is your 'searchBox' specification correct?"
                            << endl;
                    }
                }
                else
                {
                    // Test for having voxels already marked as patch
                    // -> voxel mesh is too coarse
                    if
                    (
                        voxelMeshSearch::overlaps
                        (
                            bb,
                            nDivs,
                            faceBb,
                            patchTypes,
                            static_cast<unsigned int>(patchCellType::PATCH)
                        )
                    )
                    {
                        // Determine number of voxels from number of cells
                        // in mesh
                        const labelVector& dim = mesh.geometricD();
                        forAll(dim, i)
                        {
                            if (dim[i] != -1)
                            {
                                nDivs[i] *= 2;
                            }
                        }

                        Pout<< "Voxel mesh too coarse. Bounding box "
                            << faceBb
                            << " contains both non-overset and overset patches"
                            << ". Refining voxel mesh to " << nDivs << endl;

                        return false;
                    }

                    voxelMeshSearch::fill
                    (
                        patchTypes,
                        bb,
                        nDivs,
                        faceBb,
                        patchCellType::OVERSET
                    );
                }
            }
        }
    }

    return true;
}


void Foam::cellCellStencils::trackingInverseDistance::markPatchesAsHoles
(
    PstreamBuffers& pBufs,
    const List<treeBoundBoxList>& patchBb,
    const List<labelVector>& patchDivisions,
    const PtrList<PackedList<2>>& patchParts,

    const label srcI,
    const label tgtI,
    labelList& allCellTypes
) const
{
    const pointField& allPoints = mesh_.points();
    const labelListList& allCellPoints = mesh_.cellPoints();

    const treeBoundBoxList& srcPatchBbs = patchBb[srcI];
    const treeBoundBoxList& tgtPatchBbs = patchBb[tgtI];
    const labelList& tgtCellMap = meshParts_[tgtI].cellMap();

    // 1. do processor-local src-tgt patch overlap
    {
        const treeBoundBox& srcPatchBb = srcPatchBbs[Pstream::myProcNo()];
        const treeBoundBox& tgtPatchBb = tgtPatchBbs[Pstream::myProcNo()];

        if (srcPatchBb.overlaps(tgtPatchBb))
        {
            const PackedList<2>& srcPatchTypes = patchParts[srcI];
            const labelVector& srcDivs = patchDivisions[srcI];

            forAll(tgtCellMap, tgtCelli)
            {
                label celli = tgtCellMap[tgtCelli];
                boundBox cBb(allPoints, allCellPoints[celli], false);
                cBb.min() -= smallVec_;
                cBb.max() += smallVec_;

                if
                (
                    voxelMeshSearch::overlaps
                    (
                        srcPatchBb,
                        srcDivs,
                        cBb,
                        srcPatchTypes,
                        static_cast<unsigned int>(patchCellType::PATCH)
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
                const labelVector srcDivs(is);
                const PackedList<2> srcPatchTypes(is);

                forAll(tgtCellMap, tgtCelli)
                {
                    label celli = tgtCellMap[tgtCelli];
                    boundBox cBb(allPoints, allCellPoints[celli], false);
                    cBb.min() -= smallVec_;
                    cBb.max() += smallVec_;

                    if
                    (
                        voxelMeshSearch::overlaps
                        (
                            srcPatchBb,
                            srcDivs,
                            cBb,
                            srcPatchTypes,
                            static_cast<unsigned int>(patchCellType::PATCH)
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


void Foam::cellCellStencils::trackingInverseDistance::markDonors
(
    PstreamBuffers& pBufs,
    const List<treeBoundBoxList>& meshBb,
    const PtrList<voxelMeshSearch>& meshSearches,
    const labelList& allCellTypes,

    const label srcI,
    const label tgtI,
    labelListList& allStencil,
    labelList& allDonor
) const
{
    const treeBoundBoxList& srcBbs = meshBb[srcI];
    const treeBoundBoxList& tgtBbs = meshBb[tgtI];

    const fvMesh& srcMesh = meshParts_[srcI].subMesh();
    const labelList& srcCellMap = meshParts_[srcI].cellMap();
    const voxelMeshSearch& meshSearch = meshSearches[srcI];
    const fvMesh& tgtMesh = meshParts_[tgtI].subMesh();
    const pointField& tgtCc = tgtMesh.cellCentres();
    const labelList& tgtCellMap = meshParts_[tgtI].cellMap();

    // 1. do processor-local src/tgt overlap
    {
        forAll(tgtCellMap, tgtCelli)
        {
            label srcCelli = meshSearch.findCell(tgtCc[tgtCelli]);
            if (srcCelli != -1 && allCellTypes[srcCellMap[srcCelli]] != HOLE)
            {
                label celli = tgtCellMap[tgtCelli];

                // TBD: check for multiple donors. Maybe better one? For
                //      now check 'nearer' mesh
                if (betterDonor(tgtI, allDonor[celli], srcI))
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] =
                        globalCells_.toGlobal(srcCellMap[srcCelli]);
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
            label srcCelli = meshSearch.findCell(samples[sampleI]);
            if (srcCelli != -1 && allCellTypes[srcCellMap[srcCelli]] != HOLE)
            {
                donors[sampleI] = globalCells_.toGlobal(srcCellMap[srcCelli]);
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

                // TBD: check for multiple donors. Maybe better one?
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::trackingInverseDistance::trackingInverseDistance
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    inverseDistance(mesh, dict, false),
    globalCells_(mesh_.nCells())
{
    // Initialise donor cell
    globalDonor_.setSize(mesh_.nCells());
    globalDonor_ = -1;

    // Initialise mesh partitions
    const labelIOList& zoneID = this->zoneID();
    label nZones = gMax(zoneID)+1;

    labelList nCellsPerZone(nZones, Zero);
    forAll(zoneID, celli)
    {
        nCellsPerZone[zoneID[celli]]++;
    }
    Pstream::listCombineGather(nCellsPerZone, plusEqOp<label>());
    Pstream::listCombineScatter(nCellsPerZone);

    meshParts_.setSize(nZones);
    forAll(meshParts_, zonei)
    {
        meshParts_.set
        (
            zonei,
            new fvMeshSubset(mesh_)
        );
        meshParts_[zonei].setLargeCellSubset(zoneID, zonei);

        // Trigger early evaluation of mesh dimension
        // (in case there are locally zero cells in mesh)
        (void)meshParts_[zonei].subMesh().nGeometricD();
    }


    // Print a bit
    {
        Info<< typeName << " : detected " << nZones
            << " mesh regions" << endl;
        Info<< incrIndent;
        forAll(nCellsPerZone, zonei)
        {
            Info<< indent<< "zone:" << zonei
                << " nCells:" << nCellsPerZone[zonei]
                << endl;
        }
        Info<< decrIndent;
    }


    // Do geometry update
    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::trackingInverseDistance::~trackingInverseDistance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCellStencils::trackingInverseDistance::update()
{
    DebugInfo<< FUNCTION_NAME << " : Start of analysis" << endl;

    scalar layerRelax(dict_.lookupOrDefault("layerRelax", 1.0));
    const labelIOList& zoneID = this->zoneID();
    label nZones = meshParts_.size();

    // Update stored mesh partitions for geometry changes
    forAll(meshParts_, zonei)
    {
        pointField subPoints(mesh_.points(), meshParts_[zonei].pointMap());

        fvMesh& subMesh = meshParts_[zonei].subMesh();
        subMesh.movePoints(subPoints);
    }

    DebugInfo<< FUNCTION_NAME << " : Moved zone sub-meshes" << endl;

    // Calculate fast search structure for each zone
    List<labelVector> searchBoxDivisions(dict_.lookup("searchBoxDivisions"));
    if (searchBoxDivisions.size() != nZones)
    {
        FatalIOErrorInFunction(dict_) << "Number of searchBoxDivisions "
            << searchBoxDivisions.size()
            << " should  equal the number of zones " << nZones
            << exit(FatalIOError);
    }

    const boundBox& allBb(mesh_.bounds());

    List<treeBoundBoxList> meshBb(nZones);

    // Determine zone meshes and bounding boxes
    {
        // Per processor, per zone the bounding box
        List<treeBoundBoxList> procBb(Pstream::nProcs());
        procBb[Pstream::myProcNo()].setSize(nZones);

        forAll(meshParts_, zonei)
        {
            const fvMesh& subMesh = meshParts_[zonei].subMesh();

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
        forAll(meshBb, zonei)
        {
            treeBoundBoxList& bbs = meshBb[zonei];
            bbs.setSize(Pstream::nProcs());
            forAll(procBb, proci)
            {
                bbs[proci] = procBb[proci][zonei];
            }
        }
    }

    DebugInfo<< FUNCTION_NAME << " : Calculated bounding boxes" << endl;


    // Determine patch bounding boxes. These are either global and provided
    // by the user or processor-local as a copy of the mesh bounding box.

    List<treeBoundBoxList> patchBb(nZones);
    List<labelVector> patchDivisions(searchBoxDivisions);
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

    forAll(patchParts, zonei)
    {
        while (true)
        {
            patchParts.set
            (
                zonei,
                new PackedList<2>(cmptProduct(patchDivisions[zonei]))
            );

            bool ok = markBoundaries
            (
                meshParts_[zonei].subMesh(),
                smallVec_,

                patchBb[zonei][Pstream::myProcNo()],
                patchDivisions[zonei],
                patchParts[zonei],

                meshParts_[zonei].cellMap(),
                allPatchTypes
            );

            //if (returnReduce(ok, andOp<bool>()))
            if (ok)
            {
                break;
            }
        }
    }
    DebugInfo<< FUNCTION_NAME << " : Calculated boundary voxel meshes" << endl;


    PtrList<voxelMeshSearch> meshSearches(meshParts_.size());
    forAll(meshParts_, zonei)
    {
        meshSearches.set
        (
            zonei,
            new voxelMeshSearch
            (
                meshParts_[zonei].subMesh(),
                patchBb[zonei][Pstream::myProcNo()],
                patchDivisions[zonei],
                true
            )
        );
    }

    DebugInfo<< FUNCTION_NAME << " : Constructed cell voxel meshes" << endl;

    PtrList<fvMesh> voxelMeshes;
    if (debug)
    {
        voxelMeshes.setSize(meshSearches.size());
        forAll(meshSearches, zonei)
        {
            IOobject meshIO
            (
                word("voxelMesh" + Foam::name(zonei)),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            );

            Pout<< "Constructing voxel mesh " << meshIO.objectPath() << endl;
            voxelMeshes.set(zonei, meshSearches[zonei].makeMesh(meshIO));
        }
    }


    if (debug)
    {
        forAll(patchParts, zonei)
        {
            const labelList voxelTypes(patchParts[zonei].values());
            const fvMesh& vm = voxelMeshes[zonei];
            tmp<volScalarField> tfld
            (
                createField<label>
                (
                    vm,
                    "patchTypes",
                    voxelTypes
                )
            );
            tfld().write();
        }
    }

    // Current best guess for cells
    labelList allCellTypes(mesh_.nCells(), CALCULATED);
    labelListList allStencil(mesh_.nCells());
    // zoneID of donor
    labelList allDonorID(mesh_.nCells(), -1);

    const globalIndex globalCells(mesh_.nCells());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    DebugInfo<< FUNCTION_NAME << " : Allocated donor-cell structures" << endl;

    for (label srci = 0; srci < meshParts_.size()-1; srci++)
    {
        for (label tgti = srci+1; tgti < meshParts_.size(); tgti++)
        {
            markPatchesAsHoles
            (
                pBufs,

                patchBb,
                patchDivisions,
                patchParts,

                srci,
                tgti,
                allCellTypes
            );
            markPatchesAsHoles
            (
                pBufs,

                patchBb,
                patchDivisions,
                patchParts,

                tgti,
                srci,
                allCellTypes
            );
        }
    }

    for (label srci = 0; srci < meshParts_.size()-1; srci++)
    {
        for (label tgti = srci+1; tgti < meshParts_.size(); tgti++)
        {
            markDonors
            (
                pBufs,
                meshBb,
                meshSearches,
                allCellTypes,   // to exclude hole donors

                tgti,
                srci,
                allStencil,
                allDonorID
            );
            markDonors
            (
                pBufs,
                meshBb,
                meshSearches,
                allCellTypes,   // to exclude hole donors

                srci,
                tgti,
                allStencil,
                allDonorID
            );
        }
    }

    DebugInfo<< FUNCTION_NAME << " : Determined holes and donor-acceptors"
        << endl;

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes", allCellTypes)
        );
        tfld().write();
    }


    // Use the patch types and weights to decide what to do
    forAll(allPatchTypes, celli)
    {
        if (allCellTypes[celli] != HOLE)
        {
            switch (allPatchTypes[celli])
            {
                case OVERSET:
                {
                    // Require interpolation. See if possible.

                    if (allStencil[celli].size())
                    {
                        allCellTypes[celli] = INTERPOLATED;
                    }
                    else
                    {
                        allCellTypes[celli] = HOLE;
                    }
                }
            }
        }
    }
    DebugInfo<< FUNCTION_NAME << " : Removed bad donors" << endl;

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
    findHoles(globalCells_, mesh_, zoneID, allStencil, allCellTypes);
    DebugInfo<< FUNCTION_NAME << " : Flood-filled holes" << endl;

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_hole", allCellTypes)
        );
        tfld().write();
    }

    // Add buffer interpolation layer(s) around holes
    scalarField allWeight(mesh_.nCells(), Zero);
    walkFront(layerRelax, allStencil, allCellTypes, allWeight);
    DebugInfo<< FUNCTION_NAME << " : Implemented layer relaxation" << endl;

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
    forAll(cellTypes_, celli)
    {
        if (cellTypes_[celli] == INTERPOLATED)
        {
            cellStencil_[celli].transfer(allStencil[celli]);
            cellInterpolationWeights_[celli].setSize(1);
            cellInterpolationWeights_[celli][0] = 1.0;
            interpolationCells.append(celli);
        }
        else
        {
            cellStencil_[celli].clear();
            cellInterpolationWeights_[celli].clear();
        }
    }
    interpolationCells_.transfer(interpolationCells);

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
    cellInterpolationWeight_.transfer(allWeight);
    //cellInterpolationWeight_.correctBoundaryConditions();
    dynamicOversetBlastFvMesh::correctBoundaryConditions
    <
        volScalarField,
        oversetFvPatchField<scalar>
    >(cellInterpolationWeight_.boundaryFieldRef(), false);


    if (debug & 2)
    {
        // Dump stencil
        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"injectionStencil.obj");
        Pout<< typeName << " : dumping injectionStencil to "
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

    DebugInfo<< FUNCTION_NAME << " : Transferred donor to stencil" << endl;


    // Extend stencil to get inverse distance weighted neighbours
    createStencil(globalCells);
    DebugInfo<< FUNCTION_NAME << " : Extended stencil" << endl;


    if (debug & 2)
    {
        // Dump weight
        cellInterpolationWeight_.instance() = mesh_.time().timeName();
        cellInterpolationWeight_.write();

        // Dump cell types
        volScalarField volTypes
        (
            IOobject
            (
                "cellTypes",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(volTypes.internalField(), celli)
        {
            volTypes[celli] = cellTypes_[celli];
        }
        //volTypes.correctBoundaryConditions();
        dynamicOversetBlastFvMesh::correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(volTypes.boundaryFieldRef(), false);
        volTypes.write();

        // Dump stencil
        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"stencil.obj");
        Pout<< typeName << " : dumping to " << str.name() << endl;
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
    DebugInfo<< FUNCTION_NAME << " : Finished analysis" << endl;

    // Tbd: detect if anything changed. Most likely it did!
    return true;
}


// ************************************************************************* //
