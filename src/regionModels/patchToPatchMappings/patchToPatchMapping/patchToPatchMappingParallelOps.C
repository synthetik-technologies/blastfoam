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

\*---------------------------------------------------------------------------*/

#include "patchToPatchMapping.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute>
Foam::patchToPatchMapping::constructDistributionMap
(
    const labelListList& procSendIndices
)
{
    const label nProcs = Pstream::nProcs();
    const label myProcNo = Pstream::myProcNo();

    labelListList transferSizes(nProcs);
    transferSizes[myProcNo].setSize(nProcs);
    forAll(procSendIndices, proci)
    {
        transferSizes[myProcNo][proci] = procSendIndices[proci].size();
    }
    Pstream::gatherList(transferSizes);
    Pstream::scatterList(transferSizes);

    labelListList constructMap(nProcs);
    label nSends = 0;
    forAll(constructMap, proci)
    {
        const label nTransfers = transferSizes[proci][myProcNo];
        constructMap[proci].setSize(nTransfers);
        for (label i = 0; i < nTransfers; i++)
        {
            constructMap[proci][i] = nSends++;
        }
    }

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            nSends,
            labelListList(procSendIndices),
            move(constructMap)
        )
    );
}


Foam::List<Foam::labelPair> Foam::patchToPatchMapping::distributeAddressing
(
    const mapDistribute& map
)
{
    const label nProcs = Pstream::nProcs();
    const label myProcNo = Pstream::myProcNo();

    labelListList procLocalData(nProcs);
    {
        PstreamBuffers pBuffs(Pstream::commsTypes::nonBlocking);

        // Send data
        for (label proci = 0; proci < nProcs; proci++)
        {
            const labelList& sendData = map.subMap()[proci];
            if (proci != myProcNo && sendData.size())
            {
                UOPstream(proci, pBuffs)() << sendData;
            }
        }

        pBuffs.finishedSends();

        // Set local data
        procLocalData[myProcNo] = map.subMap()[myProcNo];

        // Recieve data
        for (label proci = 0; proci < nProcs; proci++)
        {
            const labelList& recvData = map.constructMap()[proci];
            if (proci != myProcNo && recvData.size())
            {
                UIPstream(proci, pBuffs)() >> procLocalData[proci];
            }
        }
    }

    // Allocate indices
    label nLocalData = 0;
    forAll(procLocalData, proci)
    {
        nLocalData += procLocalData[proci].size();
    }
    List<labelPair> localProcData(nLocalData);


    nLocalData = 0;
    forAll(localProcData, proci)
    {
        const labelList& procData = procLocalData[proci];
        forAll(procData, i)
        {
            localProcData[nLocalData++] = {proci, procData[i]};
        }
    }

    return localProcData;
}


Foam::labelListList Foam::patchToPatchMapping::procSendIndices
(
    const labelListList& localSrcToTgtData,
    const List<labelPair>& localTgtProcData
)
{
    const label nProcs = Pstream::nProcs();
    List<labelHashSet> resSet(nProcs);

    forAll(localSrcToTgtData, datai)
    {
        const label proci = localTgtProcData[datai].first();
        resSet[proci].insert(localSrcToTgtData[datai]);
    }

    labelListList res(nProcs);
    forAll(resSet, proci)
    {
        res[proci] = resSet[proci].toc();
    }
    return res;
}

// Same as above
Foam::labelListList Foam::patchToPatchMapping::procSendIndices
(
    const List<DynamicList<label>>& localSrcToTgtData,
    const List<labelPair>& localTgtProcData
)
{
    const label nProcs = Pstream::nProcs();
    List<labelHashSet> resSet(nProcs);

    forAll(localSrcToTgtData, datai)
    {
        const label proci = localTgtProcData[datai].first();
        resSet[proci].insert(localSrcToTgtData[datai]);
    }

    labelListList res(nProcs);
    forAll(resSet, proci)
    {
        res[proci] = resSet[proci].toc();
    }
    return res;
}


void Foam::patchToPatchMapping::trimDistributionMap
(
    const boolList& oldIsUsed,
    mapDistribute& map,
    labelList& oldToNew,
    labelList& newToOld
)
{
    oldToNew.resize(oldIsUsed.size());
    newToOld.resize(count(oldIsUsed, true));

    oldToNew = -1;
    newToOld = -1;

    label newi = 0;
    forAll(oldIsUsed, oldi)
    {
        if (oldIsUsed[oldi])
        {
            oldToNew[oldi] = newi;
            newToOld[newi] = oldi;
            newi++;
        }
    }

    List<boolList> globalOldIsUsed(Pstream::nProcs());
    forAll(map.constructMap(), proci)
    {
        globalOldIsUsed[proci] =
            UIndirectList<bool>(oldIsUsed, map.constructMap()[proci]);
    }

    List<boolList> globalProcOldIsUsed(Pstream::nProcs());
    Pstream::exchange<boolList, bool>(globalOldIsUsed, globalProcOldIsUsed);

    forAll(map.subMap(), proci)
    {
        label newi = 0;
        labelList& subMap = map.subMap()[proci];
        const boolList& gProcOldIsUsed = globalProcOldIsUsed[proci];
        forAll(subMap, oldi)
        {
            if (gProcOldIsUsed[oldi])
            {
                subMap[newi++] = subMap[oldi];
            }
        }
        subMap.setSize(newi);
    }

    forAll(map.constructMap(), proci)
    {
        label newi = 0;
        labelList& constructMap = map.constructMap()[proci];
        const boolList& gOldIsUsed = globalOldIsUsed[proci];
        forAll(constructMap, oldi)
        {
            if (gOldIsUsed[oldi])
            {
                constructMap[newi++] = oldToNew[constructMap[oldi]];
            }
        }
        constructMap.setSize(newi);
    }
}


Foam::List<Foam::List<Foam::labelPair>>
Foam::patchToPatchMapping::localToRemote
(
    const labelListList& indices,
    const List<labelPair>& indexToProcIndex
)
{
    List<List<labelPair>> res(indices.size());
    if (isNull(indexToProcIndex))
    {
        const label myProcNo = Pstream::myProcNo();
        forAll(indices, datai)
        {
            const labelList& inds = indices[datai];
            res[datai].setSize(inds.size());

            forAll(inds, i)
            {
                res[datai][i] = {myProcNo, inds[i]};
            }
        }
    }
    else
    {
        forAll(indices, datai)
        {
            const labelList& inds = indices[datai];
            res[datai].setSize(inds.size());

            forAll(inds, i)
            {
                res[datai][i] = indexToProcIndex[inds[i]];
            }
        }
    }
    return res;
}

// Same as above
Foam::List<Foam::List<Foam::labelPair>>
Foam::patchToPatchMapping::localToRemote
(
    const List<DynamicList<label>>& indices,
    const List<labelPair>& indexToProcIndex
)
{
    List<List<labelPair>> res(indices.size());
    if (isNull(indexToProcIndex))
    {
        const label myProcNo = Pstream::myProcNo();
        forAll(indices, datai)
        {
            const labelList& inds = indices[datai];
            res[datai].setSize(inds.size());

            forAll(inds, i)
            {
                res[datai][i] = {myProcNo, inds[i]};
            }
        }
    }
    else
    {
        forAll(indices, datai)
        {
            const labelList& inds = indices[datai];
            res[datai].setSize(inds.size());

            forAll(inds, i)
            {
                res[datai][i] = indexToProcIndex[inds[i]];
            }
        }
    }
    return res;
}


void Foam::patchToPatchMapping::rDistributeTgtAddressing
(
    const label tgtSize,
    const mapDistribute& tgtMap,
    const List<labelPair>& localSrcProcData,
    labelListList& localSrcToTgtData
)
{
    HashTable<label, labelPair, labelPair::Hash<>> srcProcDataToLocal;
    forAll(localSrcProcData, localSrcDatai)
    {
        srcProcDataToLocal.insert
        (
            localSrcProcData[localSrcDatai],
            localSrcDatai
        );
    }

    List<List<labelPair>> srcProcToTgtData(localToRemote(localSrcToTgtData));

    rDistributeListList(tgtSize, tgtMap, srcProcToTgtData);

    localSrcToTgtData.setSize(tgtSize);
    forAll(srcProcToTgtData, tgtDatai)
    {
        const List<labelPair>& srcProcToTgt = srcProcToTgtData[tgtDatai];
        labelList& localSrcToTgt = localSrcToTgtData[tgtDatai];
        localSrcToTgt.setSize(srcProcToTgt.size());

        forAll(srcProcToTgt, i)
        {
            localSrcToTgt[i] = srcProcDataToLocal[srcProcToTgt[i]];
        }
    }
}


// Same as above
void Foam::patchToPatchMapping::rDistributeTgtAddressing
(
    const label tgtSize,
    const mapDistribute& tgtMap,
    const List<labelPair>& localSrcProcData,
    List<DynamicList<label>>& localSrcToTgtData
)
{
    HashTable<label, labelPair, labelPair::Hash<>> srcProcDataToLocal;
    forAll(localSrcProcData, localSrcDatai)
    {
        srcProcDataToLocal.insert
        (
            localSrcProcData[localSrcDatai],
            localSrcDatai
        );
    }

    List<List<labelPair>> srcProcToTgtData(localToRemote(localSrcToTgtData));

    rDistributeListList(tgtSize, tgtMap, srcProcToTgtData);

    localSrcToTgtData.setSize(tgtSize);
    forAll(srcProcToTgtData, tgtDatai)
    {
        const List<labelPair>& srcProcToTgt = srcProcToTgtData[tgtDatai];
        labelList& localSrcToTgt = localSrcToTgtData[tgtDatai];
        localSrcToTgt.setSize(srcProcToTgt.size());

        forAll(srcProcToTgt, i)
        {
            localSrcToTgt[i] = srcProcDataToLocal[srcProcToTgt[i]];
        }
    }
}


Foam::labelListList Foam::patchToPatchMapping::sendTgtPatch
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
) const
{
    List<List<treeBoundBox>> srcProcBbs(Pstream::nProcs());
    if (srcPatch.size())
    {
        srcProcBbs[Pstream::myProcNo()].setSize
        (
            1,
            makeBb(srcPatch, srcPts0, pointNormals, pointNormals0)
        );
    }

    Pstream::gatherList(srcProcBbs);
    Pstream::scatterList(srcProcBbs);

    List<DynamicList<label>> overlappingProcFaces(Pstream::nProcs());
    List<labelHashSet> overlappingProcPoints(Pstream::nProcs());
    forAll(tgtPatch, tgtFacei)
    {
        const treeBoundBox tgtBb
        (
            makeBb(tgtPatch[tgtFacei], tgtPatch.points(), tgtPts0)
        );
        forAll(srcProcBbs, proci)
        {
            const List<treeBoundBox>& procBbs = srcProcBbs[proci];
            forAll(procBbs, bbi)
            {
                if (procBbs[bbi].overlaps(tgtBb))
                {
                    overlappingProcFaces[proci].append(tgtFacei);
                    overlappingProcPoints[proci].insert(tgtPatch[tgtFacei]);
                    break;
                }
            }
        }
    }

    labelListList sendFaces(Pstream::nProcs());
    forAll(overlappingProcFaces, proci)
    {
        sendFaces[proci].transfer(overlappingProcFaces[proci]);
    }
    return sendFaces;
}


//- Create local patch by combining all valid processors
Foam::List<Foam::labelPair> Foam::patchToPatchMapping::distributePatch
(
    const mapDistribute& map,
    const primitivePatch& patch,
    const pointField& pts0,
    autoPtr<standAlonePatch>& localPatchPtr,
    autoPtr<pointField>& localPoints0Ptr
)
{
    const label nProcs = Pstream::nProcs();
    const label myProcNo = Pstream::myProcNo();
    const bool hasPoints0 = !isNull(pts0);

    List<labelList> procLocalFaceIs(nProcs);
    List<faceList> procLocalFaces(nProcs);
    List<pointField> procLocalPoints(nProcs);
    List<pointField> procLocalPoints0(nProcs);
    {
        PstreamBuffers pBuffs(Pstream::commsTypes::nonBlocking);

        // Send
        for (label proci = 0; proci < nProcs; proci++)
        {
            const labelList& sendFaceIs = map.subMap()[proci];
            if (proci != myProcNo && sendFaceIs.size())
            {
                uindirectPrimitivePatch subPatch
                (
                    UIndirectList<face>(patch, sendFaceIs),
                    patch.points()
                );
                UOPstream os(proci, pBuffs);

                if (hasPoints0)
                {
                    os  << sendFaceIs
                        << subPatch.localFaces()
                        << subPatch.localPoints()
                        << UIndirectList<point>(pts0, subPatch.meshPoints());
                }
                else
                {
                    os  << sendFaceIs
                        << subPatch.localFaces()
                        << subPatch.localPoints();
                }
            }
        }

        pBuffs.finishedSends();

        // local data
        {
            const labelList& sendFaceIs = map.subMap()[myProcNo];
            uindirectPrimitivePatch subPatch
            (
                UIndirectList<face>(patch, sendFaceIs),
                patch.points()
            );
            procLocalFaceIs[myProcNo] = sendFaceIs;
            procLocalFaces[myProcNo] = subPatch.localFaces();
            procLocalPoints[myProcNo] = subPatch.localPoints();
            if (hasPoints0)
            {
                procLocalPoints0[myProcNo] =
                    UIndirectList<point>(pts0, subPatch.meshPoints())();
            }
        }

        // Recieve
        for (label proci = 0; proci < nProcs; proci++)
        {
            if (proci != myProcNo && map.constructMap()[proci].size())
            {
                UIPstream is(proci, pBuffs);
                if (hasPoints0)
                {

                    is  >> procLocalFaceIs[proci]
                        >> procLocalFaces[proci]
                        >> procLocalPoints[proci]
                        >> procLocalPoints0[proci];
                }
                else
                {
                    is  >> procLocalFaceIs[proci]
                        >> procLocalFaces[proci]
                        >> procLocalPoints[proci];
                }
            }
        }
    }

    label nLocalPoints = 0;
    label nLocalFaces = 0;
    forAll(procLocalFaces, proci)
    {
        nLocalPoints += procLocalPoints[proci].size();
        nLocalFaces += procLocalFaces[proci].size();
    }

    List<labelPair> localProcFaces(nLocalFaces);
    faceList localFaces(nLocalFaces);
    pointField localPoints(nLocalPoints);
    pointField localPoints0(nLocalPoints);

    {
        label localPointi = 0;
        label localFacei = 0;
        forAll(procLocalFaces, proci)
        {
            const labelList& faceIs = procLocalFaceIs[proci];
            faceList& faces = procLocalFaces[proci];

            forAll(faceIs, i)
            {
                localProcFaces[localFacei] = {proci, faceIs[i]};
                face f(move(faces[i]));
                forAll(f, fpi)
                {
                    f[fpi] += localPointi;
                }
                localFaces[localFacei].transfer(f);
                localFacei++;
            }

            if (hasPoints0)
            {
                const pointField& points = procLocalPoints[proci];
                const pointField& points0 = procLocalPoints0[proci];
                forAll(points, i)
                {
                    localPoints[localPointi] = points[i];
                    localPoints0[localPointi] = points0[i];
                    localPointi++;
                }
            }
            else
            {
                const pointField& points = procLocalPoints[proci];
                forAll(points, i)
                {
                    localPoints[localPointi] = points[i];
                    localPointi++;
                }
            }
        }
    }

    localPatchPtr.reset
    (
        new standAlonePatch(move(localFaces), move(localPoints))
    );

    if (hasPoints0)
    {
        localPoints0Ptr.reset(new pointField(move(localPoints0)));
    }
    return localProcFaces;
}


// ************************************************************************* //
