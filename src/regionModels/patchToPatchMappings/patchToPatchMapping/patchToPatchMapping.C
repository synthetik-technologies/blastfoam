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
#include "treeDataPrimitivePatch.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToPatchMapping, 0);
    defineRunTimeSelectionTable(patchToPatchMapping, dictionary);
}

// * * * * * * * * * * * * * Static Members Functions  * * * * * * * * * * * //

Foam::label Foam::patchToPatchMapping::singleProcess
(
    const label sizeA,
    const label sizeB
)
{
    label procWithSize = 0;

    if (Pstream::parRun())
    {
        label hasSize = sizeA || sizeB;
        label nProcsWithSize = returnReduce(hasSize, sumOp<label>());

        if (nProcsWithSize == 0)
        {
            procWithSize = 0;
        }
        else if (nProcsWithSize == 1)
        {
            procWithSize =
                hasSize ? Pstream::myProcNo() : -1;
            reduce(procWithSize, maxOp<label>());
        }
        else
        {
            procWithSize = -1;
        }
    }
    return procWithSize;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMapping::patchToPatchMapping
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool reverse
)
:
    dict_(dict),
    srcPatch_(srcPatch),
    tgtPatch_(tgtPatch),
    reverse_(reverse),
    singleProcess_(-1),
    localSrcToTgt_(0),
    localTgtToSrc_(0),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    localSrcProcPtr_(nullptr),
    localTgtProcPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatchMapping::~patchToPatchMapping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchToPatchMapping> Foam::patchToPatchMapping::New
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool reverse
)
{
    const word type(dict.lookup("mappingType"));
    DebugInfo<< "Selecting patchToPatchMapping " << type << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown patchToPatchMapping type " << type
            << endl << endl
            << "Valid patchToPatchMapping types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<patchToPatchMapping>
    (
        cstrIter()
        (
            srcPatch,
            tgtPatch,
            dict.optionalSubDict(type + "Coeffs"),
            reverse
        )
    );
}


Foam::treeBoundBox Foam::patchToPatchMapping::makeBb
(
    const face& f,
    const pointField& pts
) const
{
    return treeBoundBox(pts, f);
}


Foam::treeBoundBox Foam::patchToPatchMapping::makeBb
(
    const face& f,
    const pointField& pts,
    const vectorField& pointNormals
) const
{
    const treeBoundBox bb(pts, f);

    const point c = bb.midpoint();
    const scalar l = bb.maxDim();

    return treeBoundBox(c - l*vector::one, c + l*vector::one);
}

Foam::treeBoundBox Foam::patchToPatchMapping::makeBb
(
    const primitivePatch& patch,
    const pointField& pts0
) const
{
    treeBoundBox bb(treeBoundBox::invertedBox);
    forAll(patch, i)
    {
        treeBoundBox bbi(makeBb(patch[i], patch.points()));
        bb = treeBoundBox
        (
            min(bb.min(), bbi.min()),
            max(bb.max(), bbi.max())
        );
    }

    if (!isNull(pts0))
    {
        forAll(patch, i)
        {
            treeBoundBox bbi(makeBb(patch[i], pts0));
            bb = treeBoundBox
            (
                min(bb.min(), bbi.min()),
                max(bb.max(), bbi.max())
            );
        }
    }
    return bb;
}

Foam::treeBoundBox Foam::patchToPatchMapping::makeBb
(
    const primitivePatch& patch,
    const pointField& pts0,
    const label facei
) const
{

    treeBoundBox bb(makeBb(patch[facei], patch.points()));
    if (!isNull(pts0))
    {
        treeBoundBox bbi(makeBb(patch[facei], pts0));
        bb = treeBoundBox
        (
            min(bb.min(), bbi.min()),
            max(bb.max(), bbi.max())
        );
    }
    return bb;
}

Foam::treeBoundBox Foam::patchToPatchMapping::makeBb
(
    const primitivePatch& patch,
    const pointField& pts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
) const
{
    treeBoundBox bb(treeBoundBox::invertedBox);
    forAll(patch, i)
    {
        treeBoundBox bbi(makeBb(patch[i], patch.points(), pointNormals));
        bb = treeBoundBox
        (
            min(bb.min(), bbi.min()),
            max(bb.max(), bbi.max())
        );
    }

    if (!isNull(pts0))
    {
        forAll(patch, i)
        {
            treeBoundBox bbi(makeBb(patch[i], pts0, pointNormals0));
            bb = treeBoundBox
            (
                min(bb.min(), bbi.min()),
                max(bb.max(), bbi.max())
            );
        }
    }
    return bb;
}


Foam::treeBoundBox Foam::patchToPatchMapping::makeBb
(
    const primitivePatch& patch,
    const pointField& pts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const label facei
) const
{

    treeBoundBox bb(makeBb(patch[facei], patch.points(), pointNormals));
    if (!isNull(pts0))
    {
        treeBoundBox bbi(makeBb(patch[facei], pts0, pointNormals0));
        bb = treeBoundBox
        (
            min(bb.min(), bbi.min()),
            max(bb.max(), bbi.max())
        );
    }
    return bb;
}


bool Foam::patchToPatchMapping::findOrIntersectFaces
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const label srcFacei,
    const label tgtFacei
)
{
    forAll(localSrcToTgt_[tgtFacei], i)
    {
        if (localSrcToTgt_[tgtFacei][i] == srcFacei)
        {
            return true;
        }
    }

    forAll(localTgtToSrc_[srcFacei], i)
    {
        if (localTgtToSrc_[srcFacei][i] == tgtFacei)
        {
            return true;
        }
    }

    return intersectFaces
    (
        srcPatch,
        srcPts0,
        tgtPatch,
        tgtPts0,
        pointNormals,
        pointNormals0,
        srcFacei,
        tgtFacei
    );
}

Foam::label Foam::patchToPatchMapping::intersectPatchQueue
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const bool isSrc,
    const DynamicList<labelPair>& queue,
    labelList& faceComplete,
    DynamicList<labelPair>& otherQueue,
    const labelList& otherFaceComplete,
    boolList& otherFaceQueued,
    boolList& otherFaceVisited
)
{
    const primitivePatch& otherPatch = isSrc ? tgtPatch : srcPatch;

    const faceList& otherLocalFaces = otherPatch.localFaces();
    const labelListList& otherPointFaces = otherPatch.pointFaces();
    const labelListList& otherFaceEdges = otherPatch.faceEdges();
    const labelListList& otherEdgeFaces = otherPatch.edgeFaces();

    DynamicList<label> otherNextFaces;
    DynamicList<label> otherCurrentFaces;
    DynamicList<label> otherVisitedFaces;

    DynamicList<label> otherQueuedFaces;
    label nFacesComplete = 0;

    forAll(queue, queuei)
    {
        const label facei = queue[queuei].first();
        const label otherFacei = queue[queuei].second();

        otherCurrentFaces.setSize(1);
        otherCurrentFaces[0] = otherFacei;
        otherVisitedFaces.setSize(1);
        otherVisitedFaces[0] = otherFacei;

        otherFaceVisited[otherFacei] = true;
        bool otherEdgeReached = false;

        while (otherCurrentFaces.size())
        {
            otherNextFaces.clear();
            forAll(otherCurrentFaces, ocfi)
            {
                const label otherFacej = otherCurrentFaces[ocfi];
                if
                (
                    findOrIntersectFaces
                    (
                        srcPatch,
                        srcPts0,
                        tgtPatch,
                        tgtPts0,
                        pointNormals,
                        pointNormals0,
                        isSrc ? facei : otherFacej,
                        isSrc ? otherFacej : facei
                    )
                )
                {
                    const face& otherFace = otherLocalFaces[otherFacej];
                    forAll(otherFace, ofpj)
                    {
                        const label otherPointj = otherFace[ofpj];
                        const labelList& oPointFaces =
                            otherPointFaces[otherPointj];
                        forAll(oPointFaces, opfj)
                        {
                            const label otherFacek = oPointFaces[opfj];
                            if (!otherFaceVisited[otherFacek])
                            {
                                otherFaceVisited[otherFacek] = true;
                                otherVisitedFaces.append(otherFacek);
                                otherQueuedFaces.append(otherFacek);
                            }
                        }

                        const label otherEdgej =
                            otherFaceEdges[otherFacej][ofpj];
                        if (otherEdgeFaces[otherEdgej].size() != 2)
                        {
                            otherEdgeReached = true;
                        }
                    }

                    if
                    (
                        otherFaceComplete[otherFacej] == 0
                     && !otherFaceQueued[otherFacej]
                    )
                    {
                        otherFaceQueued[otherFacej] = true;
                        otherQueuedFaces.append(otherFacej);
                        otherQueue.append({otherFacej, facei});
                    }
                }
            }
            otherCurrentFaces.transfer(otherNextFaces);
        }

        UIndirectList<bool>(otherFaceVisited, otherVisitedFaces) = false;

        if (faceComplete[facei] < 2)
        {
            faceComplete[facei] = otherEdgeReached ? 1 : 2;
            nFacesComplete += !otherEdgeReached;
        }
    }
    UIndirectList<bool>(otherFaceQueued, otherQueuedFaces) = false;
    return nFacesComplete;
}


void Foam::patchToPatchMapping::intersectPatches
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    if (srcPatch.empty() || tgtPatch.empty())
    {
        return;
    }

     // Build a search tree for the target patch
    typedef treeDataPrimitivePatch<primitivePatch> treeType;
    const treeBoundBox tgtTreeBox =
        treeBoundBox(tgtPatch.points(), tgtPatch.meshPoints()).extend(1e-4);
    indexedOctree<treeType> treeB
    (
        treeType
        (
            false,
            tgtPatch,
            indexedOctree<treeType>::perturbTol()
        ),
        tgtTreeBox,
        8,
        10,
        3
    );

    DebugInfo<< "Calculating patch intersections" << endl;

    label nSrcComplete = 0;
    label nTgtComplete = 0;
    labelList srcComplete(srcPatch.size(), 0);
    labelList tgtComplete(tgtPatch.size(), 0);
    boolList srcQueued(srcPatch.size(), false);
    boolList tgtQueued(tgtPatch.size(), false);
    boolList srcVisited(srcPatch.size(), false);
    boolList tgtVisited(tgtPatch.size(), false);

    DynamicList<labelPair> srcQueue, tgtQueue;

    label srcFacei = 0;
    label restarti = 0;
    while (srcFacei < srcPatch.size() && srcFacei != -1)
    {
        srcComplete[srcFacei] = 2;
        nSrcComplete++;

        const labelList tgtPatchSeeds
        (
            treeB.findBox
            (
                makeBb
                (
                    srcPatch,
                    srcPts0,
                    pointNormals,
                    pointNormals0,
                    srcFacei
                )
            )
        );

        if (tgtPatchSeeds.size())
        {
            DebugInfo
                << "Restart " << restarti
                << " from srcPatch at face "
                << srcPatch.faceCentres()[srcFacei] << endl;

            srcQueue.clear();
            tgtQueue.clear();

            forAll(tgtPatchSeeds, tgtSeedi)
            {
                const label tgtFacei = tgtPatchSeeds[tgtSeedi];
                srcQueue.append({srcFacei, tgtFacei});
                tgtQueue.append({tgtFacei, srcFacei});
            }

            label iteri = 0;
            while (true)
            {
                tgtQueue.clear();

                nSrcComplete +=
                    intersectPatchQueue
                    (
                        srcPatch,
                        srcPts0,
                        tgtPatch,
                        tgtPts0,
                        pointNormals,
                        pointNormals0,
                        true, // isSrc
                        srcQueue,
                        srcComplete,
                        tgtQueue,
                        tgtComplete,
                        tgtQueued,
                        tgtVisited
                    );

                if (!tgtQueue.size())
                {
                    break;
                }

                srcQueue.clear();

                nTgtComplete +=
                    intersectPatchQueue
                    (
                        srcPatch,
                        srcPts0,
                        tgtPatch,
                        tgtPts0,
                        pointNormals,
                        pointNormals0,
                        false,
                        tgtQueue,
                        tgtComplete,
                        srcQueue,
                        srcComplete,
                        srcQueued,
                        srcVisited
                    );

                if (!srcQueue.size())
                {
                    break;
                }

                iteri++;
            }

            DebugInfo
                << "Completed " << nSrcComplete << "/"
                << srcPatch.size() << " source faces and "
                << nTgtComplete << "/" << tgtPatch.size()
                << " target faces" << endl;
        }

        srcFacei = -1;
        forAll(srcComplete, i)
        {
            if (srcComplete[i] != 2)
            {
                srcFacei = i;
                break;
            }
        }
        restarti++;
    }
}

void Foam::patchToPatchMapping::initialise
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    localTgtToSrc_.setSize(srcPatch.size());
    forAll(localTgtToSrc_, i)
    {
        localTgtToSrc_[i].clear();
    }

    localSrcToTgt_.setSize(tgtPatch.size());
    forAll(localSrcToTgt_, i)
    {
        localSrcToTgt_[i].clear();
    }

    srcWeights_.setSize(srcPatch.size());
    forAll(srcWeights_, i)
    {
        srcWeights_[i].clear();
    }

    tgtWeights_.setSize(tgtPatch.size());
    forAll(tgtWeights_, i)
    {
        tgtWeights_[i].clear();
    }
}


Foam::labelList Foam::patchToPatchMapping::finaliseLocal
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    boolList localTgtFaceIsUsed(tgtPatch.size(), false);
    forAll(localTgtToSrc_, faceAi)
    {
        UIndirectList<bool>
        (
            localTgtFaceIsUsed,
            localTgtToSrc_[faceAi]
        ) = true;
    }


    labelList oldToNew, newToOld;
    trimDistributionMap
    (
        localTgtFaceIsUsed,
        tgtMapPtr_(),
        oldToNew,
        newToOld
    );


    forAll(localTgtToSrc_, faceAi)
    {
        labelList& tgtFaces = localTgtToSrc_[faceAi];
        forAll(tgtFaces, tgtFacei)
        {
            tgtFaces[tgtFacei] = oldToNew[tgtFaces[tgtFacei]];
        }
    }


    localSrcToTgt_ =
        List<DynamicList<label>>(localSrcToTgt_, newToOld);
    localTgtProcPtr_() =
        List<labelPair>(localTgtProcPtr_(), newToOld);

    return newToOld;
}


void Foam::patchToPatchMapping::distributeSrc
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0
)
{
    localSrcProcPtr_.reset
    (
        new List<labelPair>(distributeAddressing(srcMapPtr_()))
    );
}

void Foam::patchToPatchMapping::rDistributeTgt
(
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0
)
{
    rDistributeTgtAddressing
    (
        tgtPatch.size(),
        tgtMapPtr_(),
        localSrcProcPtr_(),
        localSrcToTgt_
    );
}


Foam::label Foam::patchToPatchMapping::finalise
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const transformer& tgtToSrc
 )
{
    label nCoupled = 0;
    forAll(localTgtToSrc_, i)
    {
        nCoupled += localTgtToSrc_[i].size();
    }
    forAll(localSrcToTgt_, i)
    {
        nCoupled += localSrcToTgt_[i].size();
    }
    return nCoupled;
}


Foam::labelList Foam::patchToPatchMapping::unmappedSrc() const
{
    DynamicList<label> unmapped(localTgtToSrc_.size());
    forAll(localTgtToSrc_, i)
    {
        if (!localTgtToSrc_[i].size())
        {
            unmapped.append(i);
        }
    }
    return unmapped;
}


Foam::labelList Foam::patchToPatchMapping::unmappedTgt() const
{
    DynamicList<label> unmapped(localSrcToTgt_.size());
    forAll(localSrcToTgt_, i)
    {
        if (!localSrcToTgt_[i].size())
        {
            unmapped.append(i);
        }
    }
    return unmapped;
}


Foam::labelList Foam::patchToPatchMapping::unmapped
(
    const primitivePatch& patch
) const
{
    if (&patch == &srcPatch_)
    {
        return unmappedSrc();
    }
    else if (&patch == &tgtPatch_)
    {
        return unmappedTgt();
    }

    FatalErrorInFunction
        << "Provided patch does not patch either patch provided for "
        << "mapping creation" << endl
        << abort(FatalError);

    return labelList();
}


void Foam::patchToPatchMapping::update
(
    const pointField& srcPts0,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const transformer& tgtToSrc
)
{

    cpuTime time;

    // Determine numbers of faces on both sides, report, and quit if either
    // side is empty
    const label srcTotalSize =
        returnReduce(srcPatch_.size(), sumOp<label>());
    const label tgtTotalSize =
        returnReduce(tgtPatch_.size(), sumOp<label>());
    if (srcTotalSize == 0 || tgtTotalSize == 0)
    {
        return;
    }

    const bool hasTgtPoints0 = !isNull(tgtPts0);

    // If a transformation is given then transform the target to the source
    tmpNrc<primitivePatch> ttgtPatchPtr(tgtPatch_);
    tmpNrc<pointField> ttgtPointsPtr(tgtPatch_.localPoints());
    tmpNrc<pointField> ttgtPoints0Ptr
    (
        hasTgtPoints0
      ? tgtPts0
      : NullObjectRef<pointField>()
    );
    if (!isNull(tgtToSrc))
    {
        ttgtPointsPtr = new pointField(tgtPatch_.localPoints());
        tgtToSrc.transformPosition
        (
            ttgtPointsPtr.ref(),
            ttgtPointsPtr.ref()
        );

        if (hasTgtPoints0)
        {
            ttgtPoints0Ptr = new pointField(tgtPts0);
            tgtToSrc.transformPosition
            (
                ttgtPoints0Ptr.ref(),
                ttgtPoints0Ptr.ref()
            );
        }

        ttgtPatchPtr =
            new primitivePatch
            (
                SubList<face>(tgtPatch_.localFaces(), tgtPatch_.size()),
                ttgtPointsPtr()
            );
    }
    const primitivePatch& ttgtPatch = ttgtPatchPtr();
    const pointField& ttgtPts0 = ttgtPoints0Ptr();

    Info<< indent << typeName << ": Calculating couplings between "
        << srcTotalSize << " source faces and " << tgtTotalSize
        << " target faces" << incrIndent << endl;

    // Determine if patches are present on multiple processors
    singleProcess_ = singleProcess(srcPatch_.size(), tgtPatch_.size());

    // Do intersection in serial or parallel as appropriate
    if (isSingleProcess())
    {
        // Initialise the workspace
        initialise
        (
            srcPatch_,
            srcPts0,
            ttgtPatch,
            ttgtPts0,
            pointNormals,
            pointNormals0
        );

        // Intersect the patches
        const treeBoundBox srcBb
        (
            makeBb(srcPatch_, srcPts0, pointNormals, pointNormals0)
        );
        const treeBoundBox tgtBb(makeBb(ttgtPatch, ttgtPts0));
        if (srcBb.overlaps(tgtBb))
        {
            intersectPatches
            (
                srcPatch_,
                srcPts0,
                ttgtPatch,
                ttgtPts0,
                pointNormals,
                pointNormals0
            );
        }
    }
    else
    {
        // Distribute the target patch so that everything is locally available
        // to the source. This is done based on bound boxes, so quite a lot of
        // faces will get distributed that ultimately are not used. These will
        // be filtered out after the intersection has been completed.
        tgtMapPtr_ =
            constructDistributionMap
            (
                sendTgtPatch
                (
                    srcPatch_,
                    srcPts0,
                    ttgtPatch,
                    ttgtPts0,
                    pointNormals,
                    pointNormals0
                )
            );

        localTgtProcPtr_.reset
        (
            new List<labelPair>
            (
                distributePatch
                (
                    tgtMapPtr_(),
                    ttgtPatch,
                    ttgtPts0,
                    localTgtPatchPtr_,
                    localTgtPoints0Ptr_
                )
            )
        );

        // Massage target patch into form that can be used by the serial
        // intersection interface
        const primitivePatch localTTgtPatch
        (
            SubList<face>
            (
                localTgtPatchPtr_(),
                localTgtPatchPtr_().size()
            ),
            localTgtPatchPtr_().points()
        );
        const pointField& localTTgtPoints0 =
            hasTgtPoints0
          ? localTgtPoints0Ptr_()
          : NullObjectRef<pointField>();

        // Initialise the workspace
        initialise
        (
            srcPatch_,
            srcPts0,
            localTTgtPatch,
            localTTgtPoints0,
            pointNormals,
            pointNormals0
        );

        // Intersect the patches
        if (localTTgtPatch.size())
        {
            intersectPatches
            (
                srcPatch_,
                srcPts0,
                localTTgtPatch,
                localTTgtPoints0,
                pointNormals,
                pointNormals0
            );
        }

        // Trim the local target patch
        finaliseLocal
        (
            srcPatch_,
            srcPts0,
            localTTgtPatch,
            localTTgtPoints0,
            pointNormals,
            pointNormals0
        );

        // Distribute the source patch
        srcMapPtr_ =
            constructDistributionMap
            (
                procSendIndices
                (
                    localSrcToTgt_,
                    localTgtProcPtr_()
                )
            );

        distributeSrc(srcPatch_, srcPts0);

        // Reverse distribute coupling data back to the target
        rDistributeTgt(tgtPatch_, tgtPts0);
    }

    // Finalise the intersection
    const label nCouples =
        finalise
        (
            srcPatch_,
            srcPts0,
            tgtPatch_,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    if (nCouples != 0)
    {
        Info<< indent
            << nCouples << " face couplings calculated in "
            << time.cpuTimeIncrement() << 's' << endl;
    }
    else
    {
        Info<< indent << "No couplings found" << endl;
    }

    Info<< decrIndent;
}


void Foam::patchToPatchMapping::update
(
    const vectorField& pointNormals,
    const transformer& tgtToSrc
)
{
    update
    (
        NullObjectRef<pointField>(),
        NullObjectRef<pointField>(),
        pointNormals,
        NullObjectRef<vectorField>(),
        tgtToSrc
    );
}

// ************************************************************************* //
