/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "burstModel.H"
#include "PatchTools.H"
#include "syncTools.H"
#include "indirectPrimitivePatch.H"
#include "uindirectPrimitivePatch.H"
#include "treeDataPrimitivePatch.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstModel, 0);
    defineRunTimeSelectionTable(burstModel, dictionary);
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::burstModel> Foam::burstModel::New
(
    const dictionary& dict,
    const bool coupled
)
{
    word modelType(dict.lookup("burstModel"));

    Info<< "Selecting burst model: " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown burst model "
            << modelType << endl << endl
            << "Valid burst models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, coupled);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::burstModel::regionize(const polyPatch& patch) const
{
    regionize(patch, patchRegionFaces_);
}


void Foam::burstModel::regionize
(
    const polyPatch& patch,
    List<labelList>& patchRegionFaces
) const
{
    const Map<label>& meshPointMap = patch.meshPointMap();
    boolList borderEdge(patch.nEdges(), false);

    labelList faceRegion;
    label nRegions = PatchTools::markZones(patch, borderEdge, faceRegion);

    Map<label> usedRegions;

    if (Pstream::parRun())
    {
        label nGlobalRegions = returnReduce(nRegions, sumOp<label>());

        Map<label> pointRegion(meshPointMap.size());
        {
            globalIndex globalRegion(nRegions);
            forAll(patch, fi)
            {
                const face& f = patch[fi];
                faceRegion[fi] = globalRegion.toGlobal(faceRegion[fi]);
                forAll(f, pi)
                {
                    pointRegion.insert(f[pi], faceRegion[fi]);
                }
            }
        }

        bool changing = true;
        while (changing)
        {
            changing = false;

            Map<label> oldPointRegion(pointRegion);

            DynamicList<label> regionMap(nGlobalRegions, -1);

            syncTools::syncPointMap
            (
                patch.boundaryMesh().mesh(),
                pointRegion,
                minEqOp<label>()
            );

            forAllConstIter(Map<label>, pointRegion, iter)
            {
                Map<label>::const_iterator oIter =
                    oldPointRegion.find(iter.key());
                if (oIter != oldPointRegion.cend())
                {
                    const label oldRegion = oIter();
                    if (oldRegion != iter())
                    {
                        label& region = regionMap(oldRegion);

                        region =
                            region < 0
                          ? iter()
                          : min(region, iter());
                        changing = true;
                    }
                }
            }

            reduce(changing, orOp<bool>());

            if (!changing)
            {
                break;
            }

            forAll(patch, fi)
            {
                const face& f = patch[fi];
                const label newRegion = regionMap[faceRegion[fi]];
                if (newRegion >= 0)
                {
                    faceRegion[fi] = newRegion;
                }
                forAll(f, pi)
                {
                    pointRegion.set(f[pi], faceRegion[fi]);
                }
            }
        }

        usedRegions.clear();
        nRegions = 0;
        forAll(faceRegion, fi)
        {
            if (usedRegions.insert(faceRegion[fi], nRegions))
            {
                nRegions++;
            }
        }

        List<labelList> globalUsedRegions(Pstream::nProcs());
        globalUsedRegions[Pstream::myProcNo()] = usedRegions.toc();
        Pstream::gatherList(globalUsedRegions);
        Pstream::scatterList(globalUsedRegions);

        usedRegions.clear();
        nRegions = 0;
        forAll(globalUsedRegions, proci)
        {
            const labelList& procUsedRegions = globalUsedRegions[proci];
            forAll(procUsedRegions, ri)
            {
                const label regioni = procUsedRegions[ri];
                if (usedRegions.insert(regioni, nRegions))
                {
                    nRegions++;
                }
            }
        }
    }
    else
    {
        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            usedRegions.insert(regioni, regioni);
        }
    }

    List<DynamicList<label>> pRegionFaces(nRegions);
    forAll(faceRegion, fi)
    {
        pRegionFaces[usedRegions[faceRegion[fi]]].append(fi);
    }

    patchRegionFaces.setSize(nRegions);
    forAll(patchRegionFaces, regioni)
    {
        patchRegionFaces[regioni].transfer(pRegionFaces[regioni]);
    }
}


void Foam::burstModel::regionize
(
    const polyPatch& patch1,
    const polyPatch& patch2
) const
{
    // Split patch1 in to regions
    regionize(patch1, patch1RegionFaces_);

    // Face to regionID for patch1
    labelList patch1FaceRegion(patch1.size());
    forAll(patch1RegionFaces_, regioni)
    {
        const labelList& regionFaces = patch1RegionFaces_[regioni];
        forAll(regionFaces, fi)
        {
            patch1FaceRegion[regionFaces[fi]] = regioni;
        }
    }

    // Calculate bounding boxes of patch1 regions
    List<treeBoundBox> patch1RegionBbs(patch1RegionFaces_.size());
    forAll(patch1RegionFaces_, regioni)
    {
        const labelList& regionFaces = patch1RegionFaces_[regioni];
        point minPt = point::max;
        point maxPt = point::min;
        forAll(regionFaces, fi)
        {
            const label facei = regionFaces[fi];
            patch1FaceRegion[facei] = regioni;

            const face& f = patch1[facei];
            forAll(f, pi)
            {
                const point& pt = patch1.points()[f[pi]];
                minPt = min(minPt, pt);
                maxPt = max(maxPt, pt);
            }
        }
        reduce(minPt, minOp<point>());
        reduce(maxPt, maxOp<point>());

        patch1RegionBbs[regioni].min() = minPt;
        patch1RegionBbs[regioni].max() = maxPt;
        patch1RegionBbs[regioni].inflate(1e-6);
    }

    // Check if a single region bounding box contains each face
    labelList patch2FaceRegion(patch2.size(), -1);
    forAll(patch2, facej)
    {
        label nContained = 0;
        label otherRegion = -1;
        const treeBoundBox bb2(patch2.points(), patch2[facej]);
        forAll(patch1RegionBbs, regioni)
        {
            if (patch1RegionBbs[regioni].contains(bb2))
            {
                nContained++;
                otherRegion = regioni;
            }
        }

        // Only one contains the face so use that region
        if (nContained == 1)
        {
            patch2FaceRegion[facej] = otherRegion;
        }
    }


    // Unmapped faces so send patches to other processors and find the
    // nearest face to take the regionID
    if (returnReduce(findIndex(patch2FaceRegion, -1) >= 0, orOp<bool>()))
    {
        const vectorField::subField& fc2 = patch2.faceCentres();

        // Determine if patches are present on multiple processors
        const bool singleProcess =
            patchToPatchTools::singleProcess
            (
                patch1.size(),
                patch2.size()
            );

        typedef treeDataPrimitivePatch<primitivePatch> DataType;
        typedef indexedOctree<DataType> TreeType;

        // Do intersection in serial or parallel as appropriate
        if (singleProcess)
        {
            TreeType tree1
            (
                DataType
                (
                    false,
                    patch1,
                    indexedOctree<TreeType>::perturbTol()
                ),
                treeBoundBox(patch1.points(), patch1.meshPoints()).extend(1e-4),
                8,
                10,
                3.0
            );

            forAll(patch2, facej)
            {
                // Not mapped yet
                if (patch2FaceRegion[facej] < 0)
                {
                    pointIndexHit info =
                        tree1.findNearest(fc2[facej], great);
                    if (info.hit())
                    {
                        patch2FaceRegion[facej] =
                            patch1FaceRegion[info.index()];
                    }
                }
            }
        }
        else
        {
            // Distribute the target patch so that everything is locally
            // available to the source. This is done based on bound boxes, so
            // quite a lot of faces will get distributed that ultimately are
            // not used. These will be filtered out after the intersection has
            // been completed.
            autoPtr<distributionMap> mapPtr =
                patchToPatchTools::constructDistributionMap
                (
                    facesToSend
                    (
                        patch2,
                        patch1 // patch to distribute
                    )
                );
            autoPtr<PrimitivePatch<faceList, pointField>> localPatch1Ptr;
            distributePatch(mapPtr(), patch1, localPatch1Ptr);

            // Massage target patch into form that can be used by the serial
            // intersection interface
            const primitivePatch localPatch1
            (
                SubList<face>(localPatch1Ptr(), localPatch1Ptr().size()),
                localPatch1Ptr().points()
            );

            // Tree for searching the local patch
            TreeType localTree1
            (
                DataType
                (
                    false,
                    localPatch1,
                    indexedOctree<TreeType>::perturbTol()
                ),
                treeBoundBox(localPatch1.points()).extend(1e-4),
                8,
                10,
                3.0
            );

            // Distribute the patch1 region
            labelList localPatch1FaceRegion(patch1FaceRegion);
            mapPtr->distribute(localPatch1FaceRegion);

            // Find the nearest region
            forAll(patch2, facei)
            {
                if (patch2FaceRegion[facei] < 0)
                {
                    pointIndexHit info =
                        localTree1.findNearest(fc2[facei], great);
                    if (info.hit())
                    {
                        patch2FaceRegion[facei] =
                            localPatch1FaceRegion[info.index()];
                    }
                }
            }
        }
    }

    // Split regions
    List<DynamicList<label>> pRegionFaces(patch1RegionFaces_.size());
    forAll(patch2FaceRegion, facej)
    {
        pRegionFaces[patch2FaceRegion[facej]].append(facej);
    }

    patch2RegionFaces_.setSize(patch1RegionFaces_.size());
    forAll(patch2RegionFaces_, regioni)
    {
        patch2RegionFaces_[regioni].transfer(pRegionFaces[regioni]);
    }
}


Foam::labelListList Foam::burstModel::facesToSend
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
) const
{
    // Get the bound boxes for the source patch. Just a single box for now.
    List<List<treeBoundBox>> srcPatchProcBbs(Pstream::nProcs());
    if (srcPatch.size())
    {
        srcPatchProcBbs[Pstream::myProcNo()] =
            List<treeBoundBox>
            (
                1,
                treeBoundBox(srcPatch.points(), srcPatch.meshPoints())
            );
    }
    else
    {
        srcPatchProcBbs[Pstream::myProcNo()] = List<treeBoundBox>();
    }

    // Distribute the boxes
    Pstream::gatherList(srcPatchProcBbs);
    Pstream::scatterList(srcPatchProcBbs);

    // Send a target face to a process if it overlaps the source patch box
    // for that process
    List<DynamicList<label>> resultDyn(Pstream::nProcs());
    forAll(tgtPatch, tgtFacei)
    {
        const treeBoundBox tgtFaceBb(tgtPatch.points(), tgtPatch[tgtFacei]);
        forAll(srcPatchProcBbs, proci)
        {
            forAll(srcPatchProcBbs[proci], bbi)
            {
                if (srcPatchProcBbs[proci][bbi].overlaps(tgtFaceBb))
                {
                    resultDyn[proci].append(tgtFacei);
                    break;
                }
            }
        }
    }

    // Transfer to non-dynamic storage
    labelListList result(Pstream::nProcs());
    forAll(result, proci)
    {
        result[proci].transfer(resultDyn[proci]);
    }

    return result;
}


Foam::List<Foam::remote> Foam::burstModel::distributePatch
(
    const distributionMap& map,
    const primitivePatch& patch,
    autoPtr<PrimitivePatch<faceList, pointField>>& localPatchPtr
) const
{
    static const label thisProci = Pstream::myProcNo();

    // Exchange per-processor data
    List<labelList> procLocalFaceis(Pstream::nProcs());
    List<faceList> procLocalFaces(Pstream::nProcs());
    List<pointField> procLocalPoints(Pstream::nProcs());
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const labelList& sendFaceis = map.subMap()[proci];

            if (proci != thisProci && sendFaceis.size())
            {
                uindirectPrimitivePatch subPatch
                (
                    UIndirectList<face>(patch, sendFaceis),
                    patch.points()
                );

                UOPstream(proci, pBufs)()
                    << sendFaceis
                    << subPatch.localFaces()
                    << subPatch.localPoints();
            }
        }

        pBufs.finishedSends();

        // Map local data
        {
            const labelList& sendFaceis = map.subMap()[thisProci];

            uindirectPrimitivePatch subPatch
            (
                UIndirectList<face>(patch, sendFaceis),
                patch.points()
            );

            procLocalFaceis[thisProci] = sendFaceis;
            procLocalFaces[thisProci] = subPatch.localFaces();
            procLocalPoints[thisProci] = subPatch.localPoints();
        }

        // Receive remote data
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != thisProci && map.constructMap()[proci].size())
            {
                UIPstream(proci, pBufs)()
                    >> procLocalFaceis[proci]
                    >> procLocalFaces[proci]
                    >> procLocalPoints[proci];
            }
        }
    }

    // Allocate
    List<remote> localProcFaces;
    faceList localFaces;
    pointField localPoints;
    {
        label nLocalFaces = 0, nLocalPoints = 0;
        forAll(procLocalFaceis, proci)
        {
            nLocalFaces += procLocalFaces[proci].size();
            nLocalPoints += procLocalPoints[proci].size();
        }
        localProcFaces.setSize(nLocalFaces);
        localFaces.setSize(nLocalFaces);
        localPoints.setSize(nLocalPoints);
    }

    // Construct the result
    label localFacei = 0, localPointi = 0;
    forAll(procLocalFaces, proci)
    {
        const labelList& fis = procLocalFaceis[proci];
        const faceList& fs = procLocalFaces[proci];
        forAll(fis, i)
        {
            localProcFaces[localFacei] = {proci, fis[i]};
            localFaces[localFacei] = face(fs[i] + localPointi);
            localFacei ++;
        }

        const pointField& ps = procLocalPoints[proci];
        forAll(ps, i)
        {
            localPoints[localPointi] = ps[i];
            localPointi ++;
        }
    }

    // Construct the local patch
    localPatchPtr.reset
    (
        new PrimitivePatch<faceList, pointField>
        (
            localFaces,
            localPoints
        )
    );

    return localProcFaces;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModel::burstModel
(
    const dictionary& dict,
    const bool coupled
)
:
    dict_(dict),
    partialBurst_(dict.lookup<bool>("partialBurst")),
    useDelta_
    (
        coupled
      ? dict.lookupOrDefault<bool>("useDelta", true)
      : false
    ),
    useAverage_(dict.lookupOrDefault<bool>("useAverage", false)),
    regionize_(dict.lookupOrDefault<bool>("regionize", !partialBurst_)),
    log_(dict.lookupOrDefault("logBurst", false)),
    needUpdate_(true),
    needRegionUpdate_(true)
{
    if (useAverage_ && partialBurst_)
    {
        WarningInFunction
            << "If partial burst is used, \"useAverage\" is ignored" << endl;
        useAverage_ = false;
    }

    if (coupled)
    {
        mappingPtr_ =
            patchToPatch::New
            (
                dict_.lookup<word>("patchToPatch"),
                false
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModel::~burstModel()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::burstModel::needUpdate() const
{
    needUpdate_ = true;
    needRegionUpdate_ = true;
}


void Foam::burstModel::updateMapping
(
    const polyPatch& patch1,
    const polyPatch& patch2
) const
{
    if (!needUpdate_ && mappingPtr_.valid())
    {
        return;
    }

    // Reset model
    mappingPtr_ = patchToPatch::New(mappingPtr_->type(), false);

    // Compute mapping
    mappingPtr_->update
    (
        patch1,
        patch1.pointNormals(),
        patch2
    );
    needUpdate_ = false;
}


bool Foam::burstModel::findBurstFaces
(
    const scalarField& pf1,
    const scalarField& magSf1,
    const scalarField& pf2,
    const scalarField& magSf2,
    const scalar burstValue,
    labelHashSet& master,
    labelHashSet& slave
) const
{
    if (useAverage_)
    {
        if (regionize_)
        {
            forAll(patch1RegionFaces_, regioni)
            {
                scalar sumPfW1 = 0.0;
                scalar sumW1 = 0.0;
                const labelList& regionFaces1 = patch1RegionFaces_[regioni];
                forAll(regionFaces1, fi)
                {
                    const label facei = regionFaces1[fi];
                    const scalar w = magSf1[facei];
                    sumPfW1 += pf1[facei]*w;
                    sumW1 += w;
                }
                const scalar pf1Mean =
                    returnReduce(sumPfW1, sumOp<scalar>())
                   /returnReduce(sumW1, sumOp<scalar>());


                scalar sumPfW2 = 0.0;
                scalar sumW2 = 0.0;
                const labelList& regionFaces2 = patch2RegionFaces_[regioni];
                forAll(regionFaces2, fi)
                {
                    const label facei = regionFaces2[fi];
                    const scalar w = magSf2[facei];
                    sumPfW2 += pf2[facei]*w;
                    sumW2 += w;
                }
                const scalar pf2Mean =
                    returnReduce(sumPfW2, sumOp<scalar>())
                   /returnReduce(sumW2, sumOp<scalar>());

                const bool regionBurst =
                    (useDelta_ && mag(pf1Mean - pf2Mean) > burstValue)
                 || (!useDelta_ && max(pf1Mean, pf2Mean) > burstValue);

                if (regionBurst)
                {
                    forAll(regionFaces1, fi)
                    {
                        master.insert(regionFaces1[fi]);
                    }

                    forAll(regionFaces2, fi)
                    {
                        slave.insert(regionFaces2[fi]);
                    }
                }
            }
        }
        else
        {
            const scalar pf1Mean = gSum(pf1*magSf1)/gSum(magSf1);
            const scalar pf2Mean = gSum(pf2*magSf2)/gSum(magSf2);
            if
            (
                (useDelta_ && mag(pf1Mean - pf2Mean) > burstValue)
            || (!useDelta_ && max(pf1Mean, pf2Mean) > burstValue)
            )
            {
                forAll(pf1, fi)
                {
                    master.insert(fi);
                }

                forAll(pf2, fi)
                {
                    slave.insert(fi);
                }
            }
        }
    }
    else
    {
        scalarField refVal1
        (
            useDelta_
          ? mag(pf1 - mappingPtr_->tgtToSrc(pf2, pf1))
          : max(mag(pf1),  mag(mappingPtr_->tgtToSrc(pf2, pf1)))
        );

        scalarField refVal2
        (
            useDelta_
          ? mag(pf2 - mappingPtr_->srcToTgt(pf1, pf2))
          : max(mag(pf2),  mag(mappingPtr_->srcToTgt(pf1, pf2)))
        );

        if (partialBurst_)
        {
            forAll(pf1, fi)
            {
                if (refVal1[fi] > burstValue)
                {
                    master.insert(fi);
                }
            }

            forAll(pf2, fi)
            {
                if (refVal2[fi] > burstValue)
                {
                    slave.insert(fi);
                }
            }
        }
        else
        {
            if (regionize_)
            {
                forAll(patch1RegionFaces_, regioni)
                {
                    bool regionBurst = false;

                    const labelList& regionFaces1 = patch1RegionFaces_[regioni];
                    const labelList& regionFaces2 = patch2RegionFaces_[regioni];

                    forAll(regionFaces1, fi)
                    {
                        const label facei = regionFaces1[fi];
                        if (refVal1[facei] > burstValue)
                        {
                            regionBurst = true;
                            break;
                        }
                    }
                    reduce(regionBurst, orOp<bool>());

                    if (regionBurst)
                    {
                        forAll(regionFaces1, fi)
                        {
                            master.insert(regionFaces1[fi]);
                        }

                        forAll(regionFaces2, fi)
                        {
                            slave.insert(regionFaces2[fi]);
                        }
                    }
                }
            }
            else
            {
                if
                (
                    returnReduce
                    (
                        gMaxMagSqr(refVal1) > burstValue
                     || gMaxMagSqr(refVal2) > burstValue,
                        orOp<bool>()
                    )
                )
                {
                    forAll(pf1, fi)
                    {
                        master.insert(fi);
                    }

                    forAll(pf2, fi)
                    {
                        slave.insert(fi);
                    }
                }
            }
        }
    }

    return
        returnReduce
        (
            master.size() + slave.size() == pf1.size() + pf2.size(),
            orOp<bool>()
        );
}


bool Foam::burstModel::findBurstFaces
(
    const scalarField& pf,
    const scalarField& magSf,
    const scalar burstValue,
    labelHashSet& faces
) const
{
    if (useAverage_)
    {
        const scalar pfMean = gSum(pf*magSf)/gSum(magSf);
        if (regionize_)
        {
            forAll(patchRegionFaces_, regioni)
            {
                const labelList& regionFaces = patchRegionFaces_[regioni];

                scalar sumPfW = 0.0;
                scalar sumW = 0.0;
                forAll(regionFaces, fi)
                {
                    const label facei = regionFaces[fi];
                    sumPfW += pf[facei]*magSf[facei];
                    sumW += magSf[facei];
                }

                const bool regionBurst =
                    returnReduce(sumPfW, sumOp<scalar>())
                   /returnReduce(sumW, sumOp<scalar>())
                  > burstValue;
                if (regionBurst)
                {
                    forAll(regionFaces, fi)
                    {
                        faces.insert(regionFaces[fi]);
                    }
                }
            }
        }
        else
        {
            if (mag(pfMean) > burstValue)
            {
                forAll(pf, facei)
                {
                    faces.insert(facei);
                }
            }
        }
    }
    else
    {
        if (partialBurst_)
        {
            forAll(pf, fi)
            {
                if (pf[fi] > burstValue)
                {
                    faces.insert(fi);
                }
            }
        }
        else
        {
            if (regionize_)
            {
                forAll(patchRegionFaces_, regioni)
                {
                    const labelList& regionFaces = patchRegionFaces_[regioni];

                    bool regionBurst = true;
                    forAll(regionFaces, fi)
                    {
                        if (pf[regionFaces[fi]] > burstValue)
                        {
                            regionBurst = true;
                            break;
                        }
                    }
                    reduce(regionBurst, orOp<bool>());

                    if (regionBurst)
                    {
                        forAll(regionFaces, fi)
                        {
                            faces.insert(regionFaces[fi]);
                        }
                    }
                }
            }
            else
            {
                if (gMax(mag(pf)) > burstValue)
                {
                    forAll(pf, facei)
                    {
                        faces.insert(facei);
                    }
                }
            }
        }
    }

    return returnReduce(faces.size() == pf.size(), andOp<bool>());
}


void Foam::burstModel::writeData(Ostream& os) const
{
    writeEntry(os, "burstModel", type());
    writeEntry(os, "partialBurst", partialBurst_);
    writeEntry(os, "useDelta", useDelta_);
    writeEntry(os, "useAverage", useAverage_);
    writeEntry(os, "regionize", regionize_);
}


// ************************************************************************* //
