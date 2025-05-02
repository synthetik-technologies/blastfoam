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


void Foam::patchToPatchMapping::sendTgtPatch
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    labelListList& sendPoints,
    labelListList& sendFaces
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

    sendFaces.setSize(Pstream::nProcs());
    sendPoints.setSize(Pstream::nProcs());
    forAll(overlappingProcFaces, proci)
    {
        sendFaces[proci].transfer(overlappingProcFaces[proci]);
        sendPoints[proci] = overlappingProcPoints[proci].toc();
    }
}


//- Create local patch by combining all valid processors
Foam::List<Foam::remote> Foam::patchToPatchMapping::distributePatch
(
    const distributionMap& map,
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

    List<remote> localProcFaces(nLocalFaces);
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


//- Create local patch by combining all valid processors
void Foam::patchToPatchMapping::distributePatch
(
    const distributionMap& map,
    const primitivePatch& patch,
    const pointField& pts0,
    List<remote>& localProcPoints,
    List<remote>& localProcFaces,
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

    localProcPoints.setSize(nLocalPoints);
    localProcFaces.setSize(nLocalFaces);
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
                    localProcPoints[localPointi] = {proci, i};
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
                    localProcPoints[localPointi] = {proci, i};
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
}

// ************************************************************************* //
