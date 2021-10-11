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

#include "box.H"
#include "mapDistribute.H"
#include "OFstream.H"
#include "meshTools.H"

const Foam::label Foam::processorLODs::box::DROP = 0;
const Foam::label Foam::processorLODs::box::REFINE = 1;
const Foam::label Foam::processorLODs::box::FIXED = 2;
const Foam::label Foam::processorLODs::box::nStartUpIter = 2;

namespace Foam
{
namespace processorLODs
{
    defineTypeNameAndDebug(box, 0);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::processorLODs::box::writeBoxes
(
    const List<DynamicList<treeBoundBox>>& fixedBoxes,
    const label iter
) const
{
    static label time = 0;

    OFstream os
    (
        "processor" + Foam::name(Pstream::myProcNo())
        + "_time" + Foam::name(time)
        + "_iter" + Foam::name(iter) + ".obj"
    );

    label verti = 0;
    for (int proci = 0; proci < Pstream::nProcs(); proci++)
//         (const int proci : Pstream::allProcs())
    {
        if (proci == Pstream::myProcNo())
        {
            continue;
        }

        const DynamicList<treeBoundBox>& procBoxes = fixedBoxes[proci];
        forAll(procBoxes, boxi)
        {
            const treeBoundBox& bb = procBoxes[boxi];

            // Write the points
            const pointField pts(bb.points());
            for (const point& p : pts)
            {
                meshTools::writeOBJ(os, p);
            }

            // Write the box faces
            for (const face& f : bb.faces)
            {
                os  << 'f';
                for (const label fpi : f)
                {
                    os  << ' ' << fpi + verti + 1;
                }
                os  << nl;
            }
            verti += pts.size();
        }
    }

    ++time;
}


void Foam::processorLODs::box::setRefineFlags
(
    const label refineIter,
    const label nTgtObjects,
    List<labelHashSet>& fixedSendElems,
    List<List<labelList>>& localTgtElems,
    List<labelList>& refineFlags,
    labelList& nElems
) const
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Identify src boxes that can be refined and send to all remote procs
    for (int proci = 0; proci < Pstream::nProcs(); proci++)
//         (const int proci : Pstream::allProcs())
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream toProc(proci, pBufs);
            toProc << nObjectsOfType_ << boxes_[proci] << newToOld_[proci];
        }
    }

    pBufs.finishedSends();

    // Test each remote src bb with local tgt objects to identify which remote
    // src boxes can/should be refined
    for (int proci = 0; proci < Pstream::nProcs(); proci++)
//         (const int proci : Pstream::allProcs())
    {
        if (proci == Pstream::myProcNo())
        {
            // Not refining boxes I send to myself - will be sending all local
            // elements
            continue;
        }

        // Receive the subset of changing src bound boxes for proci
        UIPstream fromProc(proci, pBufs);
        const label nObjects = readLabel(fromProc);
        List<treeBoundBox> remoteSrcBoxes(fromProc);
        const List<label> newToOld(fromProc);

        if (remoteSrcBoxes.empty())
        {
            continue;
        }


        labelList& procRefineFlags = refineFlags[proci];
        procRefineFlags.setSize(remoteSrcBoxes.size(), FIXED);

        if (scalar(nTgtObjects)/scalar(nObjects) < 0.1)
        {
            // Sending less than 10% of objects of this type
            // - shortcut by sending all
            fixedSendElems[proci].insert(identity(nTgtObjects));
            nElems[proci] = nTgtObjects;

            continue;
        }

        label nElem = 0;
        List<labelList>& localProcTgtElems = localTgtElems[proci];
        List<labelList> newLocalProcTgtElems(remoteSrcBoxes.size());

        forAll(remoteSrcBoxes, srcBoxi)
        {
            const treeBoundBox& remSrcBb = remoteSrcBoxes[srcBoxi];
            DynamicList<label> selectedElems(maxObjectsPerLeaf_);

            if (refineIter > nStartUpIter)
            {
                // Iterate over cached subset of tgt elements
                const label oldBoxi = newToOld[srcBoxi];
                const labelList& tgtBoxElems = localProcTgtElems[oldBoxi];

                for (const label tgtObji : tgtBoxElems)
                {
                    if (calcTgtBox(tgtObji).overlaps(remSrcBb))
                    {
                        selectedElems.append(tgtObji);
                    }
                }
            }
            else
            {
                // Iterating over all target elements
                for (label tgtObji = 0; tgtObji < nTgtObjects; ++tgtObji)
                {
                    if (calcTgtBox(tgtObji).overlaps(remSrcBb))
                    {
                        selectedElems.append(tgtObji);
                    }
                }
            }

            nElem += selectedElems.size();

            if
            (
                proci == Pstream::myProcNo()
             || selectedElems.size() < maxObjectsPerLeaf_
            )
            {
                procRefineFlags[srcBoxi] = FIXED;
                fixedSendElems[proci].insert(selectedElems);
            }
            else
            {
                procRefineFlags[srcBoxi] = REFINE;
                if (refineIter >= nStartUpIter)
                {
                    newLocalProcTgtElems[srcBoxi].transfer(selectedElems);
                }
            }
        }

        localProcTgtElems.transfer(newLocalProcTgtElems);
        nElems[proci] = nElem;
    }
}


void Foam::processorLODs::box::refineBox
(
    const label boxi,
    const label refineIter,
    const label nSrcElem,
    const treeBoundBox& origBox,
    DynamicList<treeBoundBox>& procBoxes,
    DynamicList<labelList>& procBoxElems,
    DynamicList<label>& procNewToOld
) const
{
    // Create the new boxes

    if (refineIter == nStartUpIter)
    {
        // Start caching the addressing
        for (direction octant = 0; octant < 8; ++octant)
        {
            const treeBoundBox subBb(origBox.subBbox(octant));

            // Identify the src elements for each box
            DynamicList<label> newElems(nSrcElem/2);

            for (label srcElemi = 0; srcElemi < nSrcElem; ++srcElemi)
            {
                if (subBb.overlaps(calcSrcBox(srcElemi)))
                {
                    newElems.append(srcElemi);
                }
            }

            if (newElems.size())
            {
                procBoxes.append(subBb);
                procBoxElems.append(newElems);
                procNewToOld.append(boxi);
            }
        }
    }
    else
    {
        for (direction octant = 0; octant < 8; ++octant)
        {
            const treeBoundBox subBb(origBox.subBbox(octant));

            for (label srcElemi = 0; srcElemi < nSrcElem; ++srcElemi)
            {
                if (subBb.overlaps(calcSrcBox(srcElemi)))
                {
                    procBoxes.append(subBb);
                    break;
                }
            }
        }
    }
}


void Foam::processorLODs::box::refineBox
(
    const label boxi,
    const labelList& srcAddr,
    const treeBoundBox& origBox,
    DynamicList<treeBoundBox>& procBoxes,
    DynamicList<labelList>& procBoxElems,
    DynamicList<label>& procNewToOld
) const
{
    // Create the new boxes
    for (direction octant = 0; octant < 8; ++octant)
    {
        const treeBoundBox subBb(origBox.subBbox(octant));

        // Identify the src elements for each box
        DynamicList<label> newElems(srcAddr.size()/2);

        for (const label srcElemi : srcAddr)
        {
            if (subBb.overlaps(calcSrcBox(srcElemi)))
            {
                newElems.append(srcElemi);
            }
        }

        // Only keeping new box if it overlaps src objects
        if (newElems.size())
        {
            procBoxes.append(subBb);
            procBoxElems.append(newElems);
            procNewToOld.append(boxi);
        }
    }
}


bool Foam::processorLODs::box::doRefineBoxes
(
    const label refineIter,
    const label nSrcFaces,
    const List<labelList>& refineFlags,
    const labelList& nElems,
    List<DynamicList<treeBoundBox>>& fixedBoxes
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send refine info back for list of evolving src boxes
    forAll(refineFlags, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream toProc(proci, pBufs);
            toProc << nElems[proci] << refineFlags[proci];
        }
    }

    pBufs.finishedSends();

    // Receive refine refinement actions from remote procs and use to refine
    // local src boxes
    bool refineBoxes = false;
    label nElem = 0;
    List<DynamicList<label>> newToOld(Pstream::nProcs());


    for (int proci = 0; proci < Pstream::nProcs(); proci++)
//         (const int proci : Pstream::allProcs())
    {
        if (proci == Pstream::myProcNo())
        {
            continue;
        }

        UIPstream fromProc(proci, pBufs);
        nElem += readLabel(fromProc);
        const labelList refineFlags(fromProc);

        const List<treeBoundBox>& procBoxes = boxes_[proci];
        DynamicList<treeBoundBox> newProcBoxes(procBoxes.size());
        DynamicList<labelList> newProcBoxElems(procBoxes.size());
        newToOld[proci].setCapacity(boxes_[proci].size());
        DynamicList<label>& newProcNewToOld = newToOld[proci];


        forAll(procBoxes, boxi)
        {
            if (refineFlags[boxi] == DROP)
            {
                // Skip box
            }
            else if (refineFlags[boxi] == REFINE)
            {
                if (refineIter > nStartUpIter)
                {
                    // Can use locally cached info to speed up intersection
                    // tests
                    refineBox
                    (
                        boxi,
                        boxSrcElems_[proci][boxi],
                        procBoxes[boxi],
                        newProcBoxes,
                        newProcBoxElems,
                        newProcNewToOld
                    );
                }
                else
                {
                    refineBox
                    (
                        boxi,
                        refineIter,
                        nSrcFaces,
                        procBoxes[boxi],
                        newProcBoxes,
                        newProcBoxElems,
                        newProcNewToOld
                    );
                }

                refineBoxes = true;
            }
            else if (refineFlags[boxi] == FIXED)
            {
                fixedBoxes[proci].append(procBoxes[boxi]);
            }
            else
            {
                FatalErrorInFunction
                    << "Unhandled refine action " << refineFlags[boxi]
                    << abort(FatalError);
            }
        }

        // Only keeping boxes that are to be refined
        boxes_[proci].transfer(newProcBoxes);
        boxSrcElems_[proci].transfer(newProcBoxElems);
        newToOld_[proci].transfer(newProcNewToOld);
    }

    return returnReduce(refineBoxes, orOp<bool>());
}


Foam::autoPtr<Foam::mapDistribute> Foam::processorLODs::box::createMap
(
    const label nSrcElems,
    const label nTgtElems
)
{
    // Store elements to send - will be used to build the mapDistribute
    List<labelHashSet> fixedSendElems(Pstream::nProcs());

    // List of local tgt elems to optimise searching for tgt elements inside
    // remote src boxes
    List<List<labelList>> localTgtElems(Pstream::nProcs());

    // Storage of boxes per processor - only useful for debugging
    List<DynamicList<treeBoundBox>> fixedBoxes(Pstream::nProcs());

    // Iterate to subdivide src boxes
    label refinementIter = 1;
    bool refineSrcBoxes = true;
    while (refineSrcBoxes && (refinementIter <= nRefineIterMax_))
    {
        // Per processor refinement info
        List<labelList> refineFlags(Pstream::nProcs());
        labelList nElems(Pstream::nProcs(), Zero);

        // Assess how many target elements intersect the source bounding boxes
        // and use the info to flag how the source boxes should be refined
        setRefineFlags
        (
            refinementIter,
            nTgtElems,
            fixedSendElems,
            localTgtElems,
            refineFlags,
            nElems
        );

        // Refine the source bounding boxes
        refineSrcBoxes =
            doRefineBoxes
            (
                refinementIter,
                nSrcElems,
                refineFlags,
                nElems,
                fixedBoxes
            );

        ++refinementIter;

        if (debug > 1)
        {
            // Include any boxes that are still evolving
            List<DynamicList<treeBoundBox>> allBoxes(fixedBoxes);
            forAll(allBoxes, proci)
            {
                allBoxes[proci].append(boxes_[proci]);
            }
            writeBoxes(allBoxes, refinementIter);
        }
    }

    if (debug)
    {
        Pout<< "Local src boxes after " << refinementIter-1 << " iterations:"
            << nl;

        forAll(fixedBoxes, proci)
        {
            // Include any boxes that are still evolving in box count
            label nBox = fixedBoxes[proci].size() + boxes_[proci].size();
            Pout<< "    proc:" << proci << " nBoxes:" << nBox << nl;
        }
        Pout<< endl;
    }

    // Convert send elems into a List<labelList>
    List<labelList> sendElems(Pstream::nProcs());
    forAll(localTgtElems, proci)
    {
        if (proci == Pstream::myProcNo() && nSrcElems)
        {
            sendElems[proci] = identity(nTgtElems);
        }
        else
        {
            labelHashSet& allIDs = fixedSendElems[proci];

            const List<labelList>& procBoxElems = localTgtElems[proci];

            for (const labelList& elems: procBoxElems)
            {
                allIDs.insert(elems);
            }

            sendElems[proci] = allIDs.toc();
        }
    }

    fixedSendElems.clear();
    localTgtElems.clear();

    if (debug)
    {
        Pout<< "Local objects: " << nTgtElems << " Send map:" << nl
            << tab << "proc" << tab << "objects" << endl;

        forAll(sendElems, proci)
        {
            Pout<< tab << proci << tab << sendElems[proci].size()
                << endl;
        }
    }

    return createLODMap(sendElems);
}


Foam::autoPtr<Foam::mapDistribute> Foam::processorLODs::box::createLODMap
(
    List<labelList>& sendElems
) const
{
    // Send over how many objects I need to receive
    const label localProci = Pstream::myProcNo();
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[localProci].setSize(Pstream::nProcs());
    forAll(sendElems, proci)
    {
        sendSizes[localProci][proci] = sendElems[proci].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // My local segment first
    constructMap[localProci] = identity(sendElems[localProci].size());

    label segmenti = constructMap[localProci].size();
    forAll(constructMap, proci)
    {
        if (proci != localProci)
        {
            // What I need to receive is what other processor is sending to me
            label nRecv = sendSizes[proci][localProci];
            constructMap[proci].setSize(nRecv);

            for (label& addr : constructMap[proci])
            {
                addr = segmenti++;
            }
        }
    }

    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            segmenti,                   // size after construction
            std::move(sendElems),
            std::move(constructMap)
        )
    );

    return mapPtr;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::processorLODs::box::box
(
    const UList<point>& srcPoints,
    const UList<point>& tgtPoints,
    const label maxObjectsPerLeaf,
    const label nObjectsOfType,
    const label nRefineIterMax
)
:
    processorLOD(maxObjectsPerLeaf, nObjectsOfType),
    srcPoints_(srcPoints),
    tgtPoints_(tgtPoints),
    boxes_(Pstream::nProcs()),
    nRefineIterMax_(nRefineIterMax),
    newToOld_(Pstream::nProcs()),
    boxSrcElems_(Pstream::nProcs())
{
    // Initialise each processor with a single box large enough to include all
    // of its local src points
    if (srcPoints_.size())
    {
        forAll(boxes_, proci)
        {
            List<treeBoundBox>& procBoxes = boxes_[proci];

            // Note: inflate to ease overlap operations and to handle 2-D
            // axis-aligned objects
            treeBoundBox srcBb(srcPoints_);
            srcBb.inflate(0.01);

            DynamicList<treeBoundBox> newProcBoxes(1);
            newProcBoxes.append(srcBb);
            procBoxes.transfer(newProcBoxes);
        }
    }
}


// ************************************************************************* //
