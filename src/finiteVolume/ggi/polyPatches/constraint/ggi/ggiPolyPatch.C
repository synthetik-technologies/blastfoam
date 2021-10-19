/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "ggiPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "SubField.H"
#include "Time.H"
#include "indirectPrimitivePatch.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        ggiPolyPatch,
        0
    );

    addToRunTimeSelectionTable(polyPatch, ggiPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, ggiPolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ggiPolyPatch::active() const
{
    polyPatchID shadow(shadowName_, boundaryMesh());
    faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

    if (shadow.active() && zone.active())
    {
        if (!Pstream::parRun() && !localParallel())
        {
            // Patch is present in serial run, but zone is not the same size
            // Probably doing decomposition and reconstruction
            // HJ, 14/Sep/2016
            return false;
        }

        return true;
    }
    else
    {
        // Zones not active.  Nothing to check
        return false;
    }
}


void Foam::ggiPolyPatch::calcZoneAddressing() const
{
    // Calculate patch-to-zone addressing
    if (zoneAddressingPtr_)
    {
        FatalErrorInFunction
            << "Patch to zone addressing already calculated"
            << abort(FatalError);
    }

    // Calculate patch-to-zone addressing
    zoneAddressingPtr_ = new labelList(size());
    labelList& zAddr = *zoneAddressingPtr_;
    const faceZone& myZone = zone();

    for (label i = 0; i < size(); i++)
    {
        zAddr[i] = myZone.whichFace(start() + i);
    }

    // Check zone addressing
    if (zAddr.size() > 0 && min(zAddr) < 0)
    {
        FatalErrorInFunction
            << "Problem with patch-to-zone addressing: some patch faces "
            << "not found in interpolation zone"
            << abort(FatalError);
    }
}


void Foam::ggiPolyPatch::calcRemoteZoneAddressing() const
{
    // Calculate patch-to-remote zone addressing
    if (remoteZoneAddressingPtr_)
    {
        FatalErrorInFunction
            << "Patch to remote zone addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "ggiPolyPatch::calcRemoteZoneAddressing() const for patch "
            << name() << endl;
    }

    // Once zone addressing is established, visit the opposite side and find
    // out which face data is needed for interpolation from the remote zone
    // to interpolate to my live faces
    boolList usedShadows(shadow().zone().size(), false);

    const labelList& zAddr = zoneAddressing();

    if (owner())
    {
        const labelListList& addr = patchToPatch().masterAddr();

        forAll (zAddr, mfI)
        {
            const labelList& nbrs = addr[zAddr[mfI]];

            forAll (nbrs, nbrI)
            {
                usedShadows[nbrs[nbrI]] = true;
            }
        }
    }
    else
    {
        const labelListList& addr = patchToPatch().slaveAddr();

        forAll (zAddr, mfI)
        {
            const labelList& nbrs = addr[zAddr[mfI]];

            forAll (nbrs, nbrI)
            {
                usedShadows[nbrs[nbrI]] = true;
            }
        }
    }

    // Count and pick up shadow indices
    label nShadows = 0;

    forAll (usedShadows, sI)
    {
        if (usedShadows[sI])
        {
            nShadows++;
        }
    }

    remoteZoneAddressingPtr_ = new labelList(nShadows);
    labelList& rza = *remoteZoneAddressingPtr_;

    // Reset counter for re-use
    nShadows = 0;

    forAll (usedShadows, sI)
    {
        if (usedShadows[sI])
        {
            rza[nShadows] = sI;
            nShadows++;
        }
    }
}


void Foam::ggiPolyPatch::calcPatchToPatch() const
{
    // Create patch-to-patch interpolation
    if (patchToPatchPtr_)
    {
        FatalErrorInFunction
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    if (owner())
    {
        if (debug)
        {
            InfoInFunction
                << "Calculating patch to patch interpolation for patch"
                << name() << endl;
        }

        // Create interpolation for zones
        patchToPatchPtr_ =
            new ggiZoneInterpolation
            (
                zone()(),           // This zone reference
                shadow().zone()(),  // Shadow zone reference
                forwardT(),
                reverseT(),
                -separation(), // Slave-to-master separation: Use - localValue
                true,          // Patch data is complete on all processors
                // Bug fix, delayed slave evaluation causes error
                // HJ, 30/Jun/2013
                SMALL,         // Non-overlapping face tolerances
                SMALL,         // HJ, 24/Oct/2008
                true,          // Rescale weighting factors.  Bug fix, MB.
                reject_        // Quick rejection algorithm, default BB_OCTREE
            );

        // Abort immediately if uncovered faces are present and the option
        // bridgeOverlap is not set.
        if
        (
            (
                patchToPatch().uncoveredMasterFaces().size() > 0
            && !bridgeOverlap()
            )
         || (
                patchToPatch().uncoveredSlaveFaces().size() > 0
            && !shadow().bridgeOverlap()
            )
        )
        {
            FatalErrorInFunction
                << "Found uncovered faces for GGI interface "
                << name() << "/" << shadowName()
                << " while the bridgeOverlap option is not set "
                << "in the boundary file." << endl
                << "This is an unrecoverable error. Aborting."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("void ggiPolyPatch::calcPatchToPatch() const")
            << "Attempting to create GGIInterpolation on a shadow"
            << abort(FatalError);
    }
}


void Foam::ggiPolyPatch::calcReconFaceCellCentres() const
{
    if (reconFaceCellCentresPtr_)
    {
        FatalErrorInFunction
            << "Reconstructed cell centres already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoInFunction
            << "Calculating recon centres for patch "
            << name() << endl;
    }

    // Create neighbouring face centres using interpolation
    if (owner())
    {
        const label shadowID = shadowIndex();

        // Get interpolated shadow face centre to face cell centre vectors
        tmp<vectorField> tdf = interpolate
        (
            boundaryMesh()[shadowID].faceCellCentres()
          - boundaryMesh()[shadowID].faceCentres()
        );

        // Get face centres on master side
        const vectorField::subField cf(faceCentres());

        if (bridgeOverlap_)
        {
            // Get face cell centres on master side
            const vectorField ccf(faceCellCentres());

            // Deltas for fully uncovered faces
            const vectorField uncoveredDeltas(cf - ccf);

            // Set uncovered deltas to fully uncovered faces
            setUncoveredFaces(uncoveredDeltas, tdf.ref());

            // Scale partially overlapping faces
            scalePartialFaces(tdf.ref());
        }

        // Calculate the reconstructed cell centres
        reconFaceCellCentresPtr_ = new vectorField(tdf() + cf);
    }
    else
    {
        FatalErrorInFunction
            << "Attempting to create reconFaceCellCentres on a shadow"
            << abort(FatalError);
    }
}


void Foam::ggiPolyPatch::calcLocalParallel() const
{
    // Calculate patch-to-zone addressing
    if (localParallelPtr_)
    {
        FatalErrorInFunction
            << "Local parallel switch already calculated"
            << abort(FatalError);
    }

    localParallelPtr_ = new bool(false);
    bool& emptyOrComplete = *localParallelPtr_;

    if (size() > zone().size())
    {
        FatalErrorInFunction
            << "Patch size is greater than zone size for GGI patch "
            << name() << ".  This is not allowed: "
            << "the face zone must contain all patch faces and be "
            << "global in parallel runs"
            << abort(FatalError);
    }

    // Calculate localisation on master and shadow
    if ((size() == 0 && shadow().size() == 0))
    {
        // No ggi on this processor
        emptyOrComplete = true;
    }
    else if (!zone().empty() || !shadow().zone().empty())
    {
        // GGI present on the processor and complete for both
        emptyOrComplete =
        (
            zone().size() == size()
         && shadow().zone().size() == shadow().size()
        );
    }
    else
    {
        // Master and shadow on different processors
        emptyOrComplete = false;
    }

    reduce(emptyOrComplete, andOp<bool>());

    // Note: only master allocates the comm_
    // HJ, 20/Sep/2016
    if (!emptyOrComplete && Pstream::parRun() && owner())
    {
        // Count how many patch faces exist on each processor
        labelList nFacesPerProc(Pstream::nProcs(), 0);
        nFacesPerProc[Pstream::myProcNo()] = size() + shadow().size();

        Pstream::gatherList(nFacesPerProc);
        Pstream::scatterList(nFacesPerProc);

        // Make a comm, from all processors that contain the ggi faces
        labelList ggiCommProcs(Pstream::nProcs());
        label nGgiCommProcs = 0;

        forAll (nFacesPerProc, procI)
        {
            if (nFacesPerProc[procI] > 0)
            {
                ggiCommProcs[nGgiCommProcs] = procI;
                nGgiCommProcs++;
            }
        }
        ggiCommProcs.setSize(nGgiCommProcs);

        // Allocate communicator
        comm_ = Pstream::allocateCommunicator
        (
            Pstream::worldComm,
            ggiCommProcs
        );

        if (debug && Pstream::parRun())
        {
            Info<< "Allocating communicator for GGI patch " << name()
                << " with " << ggiCommProcs << ": " << comm_
                << endl;
        }
    }
    else
    {
        comm_ = boundaryMesh().mesh().comm();
    }

    if (debug && Pstream::parRun())
    {
        Info<< "GGI patch Master: " << name()
            << " Slave: " << shadowName() << " is ";

        if (emptyOrComplete)
        {
           Info<< "local parallel" << endl;
        }
        else
        {
            Info<< "split between multiple processors" << endl;
        }
    }
}


void Foam::ggiPolyPatch::calcSendReceive() const
{
    // Note: all processors will execute calcSendReceive but only master will
    // hold complete the information.  Therefore, pointers on slave processors
    // will remain meaningless, but for purposes of consistency
    // (of the calc-call) they will be set to zero-sized array
    // HJ, 4/Jun/2011

    if (mapPtr_)
    {
        FatalErrorInFunction
            << "Send-receive addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "ggiPolyPatch::calcSendReceive() const for patch "
            << name() << endl;
    }

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // Gather send and receive addressing (to master)

    // Get patch-to-zone addressing
    const labelList& za = zoneAddressing();

    // Make a zone-sized field and fill it in with proc markings for processor
    // that holds and requires the data
    labelField zoneProcID(zone().size(), -1);

    forAll (za, zaI)
    {
        zoneProcID[za[zaI]] = Pstream::myProcNo();
    }

    reduce(zoneProcID, maxOp<labelList>());

    const labelList& shadowRza = shadow().remoteZoneAddressing();

    // Find out where my zone data is coming from
    labelList nRecv(Pstream::nProcs(), 0);

    // Note: only visit the data from the local zone
    forAll (shadowRza, shadowRzaI)
    {
        nRecv[zoneProcID[shadowRza[shadowRzaI]]]++;
    }

    // Make a receiving sub-map
    // It tells me which data I will receive from which processor and
    // where I need to put it into the remoteZone data before the mapping
    labelListList constructMap(Pstream::nProcs());

    // Size the receiving list
    forAll (nRecv, procI)
    {
        constructMap[procI].setSize(nRecv[procI]);
    }

    // Reset counters for processors
    nRecv = 0;

    forAll (shadowRza, shadowRzaI)
    {
        label recvProc = zoneProcID[shadowRza[shadowRzaI]];

        constructMap[recvProc][nRecv[recvProc]] = shadowRza[shadowRzaI];

        nRecv[recvProc]++;
    }

    // Make the sending sub-map
    // It tells me which data is required from me to be sent to which
    // processor

    // Algorithm
    // - expand the local zone faces with indices into a size of local zone
    // - go through remote zone addressing on all processors
    // - find out who hits my faces
    labelList localZoneIndices(zone().size(), -1);

    forAll (za, zaI)
    {
        localZoneIndices[za[zaI]] = zaI;
    }

    labelListList shadowToReceiveAddr(Pstream::nProcs());

    // Get the list of what my shadow needs to receive from my zone
    // on all other processors
    shadowToReceiveAddr[Pstream::myProcNo()] = shadowRza;
    Pstream::gatherList(shadowToReceiveAddr);
    Pstream::scatterList(shadowToReceiveAddr);

    // Now local zone indices contain the index of a local face that will
    // provide the data.  For faces that are not local, the index will be -1

    // Find out where my zone data is going to

    // Make a sending sub-map
    // It tells me which data I will send to which processor
    labelListList sendMap(Pstream::nProcs());

    // Collect local labels to be sent to each processor
    forAll (shadowToReceiveAddr, procI)
    {
        const labelList& curProcSend = shadowToReceiveAddr[procI];

        // Find out how much of my data is going to this processor
        label nProcSend = 0;

        forAll (curProcSend, sendI)
        {
            if (localZoneIndices[curProcSend[sendI]] > -1)
            {
                nProcSend++;
            }
        }

        if (nProcSend > 0)
        {
            // Collect the indices
            labelList& curSendMap = sendMap[procI];

            curSendMap.setSize(nProcSend);

            // Reset counter
            nProcSend = 0;

            forAll (curProcSend, sendI)
            {
                if (localZoneIndices[curProcSend[sendI]] > -1)
                {
                    curSendMap[nProcSend] =
                        localZoneIndices[curProcSend[sendI]];
                    nProcSend++;
                }
            }
        }
    }

    // Map will return the object of the size of remote zone
    // HJ, 9/May/2016
    mapPtr_ = new mapDistribute(zone().size(), move(sendMap), move(constructMap));
}


void Foam::ggiPolyPatch::clearGeom()
{
    deleteDemandDrivenData(reconFaceCellCentresPtr_);

    // Remote addressing and send-receive maps depend on the local
    // position.  Therefore, it needs to be recalculated at mesh motion.
    // Local zone addressing does not change with mesh motion
    // HJ, 23/Jun/2011
    deleteDemandDrivenData(remoteZoneAddressingPtr_);

    // localParallel depends on geometry - must be cleared!
    // HR, 11/Jul/2013
    deleteDemandDrivenData(localParallelPtr_);

    // Clear communicator
    if (comm_ != Pstream::worldComm)
    {
        Pstream::freeCommunicator(comm_);
        comm_ = Pstream::worldComm;
    }

    deleteDemandDrivenData(mapPtr_);
}


void Foam::ggiPolyPatch::clearOut()
{
    clearGeom();

    shadowIndex_ = -1;
    zoneIndex_ = -1;

    deleteDemandDrivenData(zoneAddressingPtr_);
    deleteDemandDrivenData(patchToPatchPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiPolyPatch::ggiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    seperatedPolyPatch(name, size, start, index, bm, patchType),
    shadowName_("initializeMe"),
    zoneName_("initializeMe"),
    bridgeOverlap_(false),
    reject_(ggiZoneInterpolation::BB_OCTREE),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    remoteZoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    comm_(Pstream::worldComm),
    tag_(Pstream::msgType()),
    mapPtr_(nullptr)
{}


Foam::ggiPolyPatch::ggiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& shadowName,
    const word& zoneName,
    const bool bridgeOverlap,
    const ggiZoneInterpolation::quickReject reject
)
:
    seperatedPolyPatch(name, size, start, index, bm, patchType),
    shadowName_(shadowName),
    zoneName_(zoneName),
    bridgeOverlap_(bridgeOverlap),
    reject_(reject),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    remoteZoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    comm_(Pstream::worldComm),
    tag_(Pstream::msgType()),
    mapPtr_(nullptr)
{}


Foam::ggiPolyPatch::ggiPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    seperatedPolyPatch(name, dict, index, bm, patchType),
    shadowName_(dict.lookup("shadowPatch")),
    zoneName_(dict.lookup("zone")),
    bridgeOverlap_(dict.lookup("bridgeOverlap")),
    reject_(ggiZoneInterpolation::BB_OCTREE),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    remoteZoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    comm_(Pstream::worldComm),
    tag_(Pstream::msgType()),
    mapPtr_(nullptr)
{
    if (dict.found("quickReject"))
    {
        reject_ = ggiZoneInterpolation::quickRejectNames_.read
        (
            dict.lookup("quickReject")
        );
    }
}


// Construct as copy, resetting the face list and boundary mesh data
Foam::ggiPolyPatch::ggiPolyPatch
(
    const ggiPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    seperatedPolyPatch(pp, bm, index, newSize, newStart),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    bridgeOverlap_(pp.bridgeOverlap_),
    reject_(pp.reject_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    remoteZoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    comm_(pp.comm_),
    tag_(pp.tag_),
    mapPtr_(nullptr)
{}


Foam::ggiPolyPatch::ggiPolyPatch
(
    const ggiPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    seperatedPolyPatch(pp, bm),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    bridgeOverlap_(pp.bridgeOverlap_),
    reject_(pp.reject_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    remoteZoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    comm_(pp.comm_),
    tag_(pp.tag_),
    mapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ggiPolyPatch::~ggiPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::ggiPolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1 && shadowName_ != Foam::word::null)
    {
        // Grab shadow patch index
        polyPatchID shadow(shadowName_, boundaryMesh());

        if (!shadow.active())
        {
            FatalErrorIn("label ggiPolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadowName_
                << " not found.  Please check your GGI interface definition."
                << abort(FatalError);
        }

        shadowIndex_ = shadow.index();

        // Check the other side is a ggi
        if (!isA<ggiPolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label ggiPolyPatch::shadowIndex() const")
                << "Shadow of ggi patch " << name()
                << " named " << shadowName() << " is not a ggi.  Type: "
                << boundaryMesh()[shadowIndex_].type() << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }

        // Check for GGI onto self
        if (index() == shadowIndex_)
        {
            FatalErrorIn("label ggiPolyPatch::shadowIndex() const")
                << "ggi patch " << name() << " created as its own shadow"
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}


Foam::label Foam::ggiPolyPatch::zoneIndex() const
{
    if (zoneIndex_ == -1 && zoneName_ != Foam::word::null)
    {
        // Grab zone patch index
        faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

        if (!zone.active())
        {
            FatalErrorIn("label ggiPolyPatch::zoneIndex() const")
                << "Face zone name " << zoneName_
                << " for GGI patch " << name()
                << " not found.  Please check your GGI interface definition."
                << abort(FatalError);
        }

        zoneIndex_ = zone.index();
    }

    return zoneIndex_;
}


const Foam::ggiPolyPatch& Foam::ggiPolyPatch::shadow() const
{
    return refCast<const ggiPolyPatch>(boundaryMesh()[shadowIndex()]);
}


const Foam::faceZone& Foam::ggiPolyPatch::zone() const
{
    return boundaryMesh().mesh().faceZones()[zoneIndex()];
}


Foam::label Foam::ggiPolyPatch::comm() const
{
    // Note: comm is calculated with localParallel and will use the
    // localParallelPtr_ for signalling.  HJ, 10/Sep/2016
    if (owner())
    {
        if (!localParallelPtr_)
        {
            calcLocalParallel();
        }

        return comm_;
    }
    else
    {
        return shadow().comm();
    }
}


int Foam::ggiPolyPatch::tag() const
{
    return Pstream::msgType();
}


const Foam::labelList& Foam::ggiPolyPatch::zoneAddressing() const
{
    if (!zoneAddressingPtr_)
    {
        calcZoneAddressing();
    }

    return *zoneAddressingPtr_;
}


const Foam::labelList& Foam::ggiPolyPatch::remoteZoneAddressing() const
{
    if (!remoteZoneAddressingPtr_)
    {
        calcRemoteZoneAddressing();
    }

    return *remoteZoneAddressingPtr_;
}


bool Foam::ggiPolyPatch::localParallel() const
{
    // Calculate patch-to-zone addressing
    if (!localParallelPtr_)
    {
        calcLocalParallel();
    }

    return *localParallelPtr_;
}


const Foam::ggiZoneInterpolation& Foam::ggiPolyPatch::patchToPatch() const
{
    if (owner())
    {
        if (!patchToPatchPtr_)
        {
            Info<< "Initializing the GGI interpolator between "
                << "master/shadow patches: "
                << name() << "/" << shadowName()
                << endl;

            calcPatchToPatch();
        }

        return *patchToPatchPtr_;
    }
    else
    {
        return shadow().patchToPatch();
    }
}


const Foam::mapDistribute& Foam::ggiPolyPatch::map() const
{
    if (!mapPtr_)
    {
        calcSendReceive();
    }

    return *mapPtr_;
}


const Foam::vectorField& Foam::ggiPolyPatch::reconFaceCellCentres() const
{
    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }

    return *reconFaceCellCentresPtr_;
}


// void Foam::ggiPolyPatch::initAddressing()
// {
//     if (active())
//     {
//         // Calculate transforms for correct GGI cut
//         calcTransforms();
//
//         if (owner())
//         {
//             shadow().calcTransforms();
//         }
//
//         // Force zone addressing and remote zone addressing
//         // (uses GGI interpolator)
//         zoneAddressing();
//         remoteZoneAddressing();
//
//         // Force local parallel
//         if (Pstream::parRun() && !localParallel())
//         {
//             // Calculate send addressing
//             map();
//         }
//     }
//
//     polyPatch::initAddressing(buf);
// }


// void Foam::ggiPolyPatch::calcAddressing()
// {
//     polyPatch::calcAddressing();
// }


void Foam::ggiPolyPatch::initCalcGeometry(PstreamBuffers& buf)
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011
    if (active())
    {
        // Note: Only master calculates recon; slave uses master interpolation
        if (owner())
        {
            reconFaceCellCentres();
        }
    }

    polyPatch::initCalcGeometry(buf);
}


void Foam::ggiPolyPatch::calcGeometry(PstreamBuffers& buf)
{
    polyPatch::calcGeometry(buf);

    // Note: Calculation of transforms must be forced before the
    // reconFaceCellCentres in order to correctly set the transformation
    // in the interpolation routines.
    // HJ, 3/Jul/2009
}


void Foam::ggiPolyPatch::initMovePoints(PstreamBuffers& buf, const pointField& p)
{
    clearGeom();

    // Calculate transforms on mesh motion?
    calcTransforms();

    if (owner())
    {
        const_cast<ggiPolyPatch&>(shadow()).clearGeom();
        shadow().calcTransforms();
    }

    // Update interpolation for new relative position of GGI interfaces
    if (patchToPatchPtr_)
    {
        patchToPatchPtr_->movePoints
        (
            forwardT(),
            reverseT(),
            -separation()
        );
    }

    // Recalculate send and receive maps
    if (active())
    {
        // Force zone addressing first
        zoneAddressing();
        remoteZoneAddressing();

        if (Pstream::parRun() && !localParallel())
        {
            map();
        }

        if (owner())
        {
            reconFaceCellCentres();
        }
    }

    polyPatch::initMovePoints(buf, p);
}


void Foam::ggiPolyPatch::movePoints(PstreamBuffers& buf, const pointField& p)
{
    polyPatch::movePoints(buf, p);
}


void Foam::ggiPolyPatch::initUpdateMesh(PstreamBuffers& buf)
{
    polyPatch::initUpdateMesh(buf);
}


void Foam::ggiPolyPatch::updateMesh(PstreamBuffers& buf)
{
    polyPatch::updateMesh(buf);
    clearOut();
}


void Foam::ggiPolyPatch::calcTransforms() const
{
    // Simplest interface: no transform or separation.  HJ, 24/Oct/2008
    forwardT_.setSize(0);
    reverseT_.setSize(0);
    separation_.setSize(0);

    if (debug > 1 && owner())
    {
//         if
//         (
//             !empty()
//          && patchToPatch().uncoveredMasterFaces().size() > 0
//         )
//         {
//             // Write uncovered master faces
//             Info<< "Writing uncovered master faces for patch "
//                 << name() << " as VTK." << endl;
//
//             const polyMesh& mesh = boundaryMesh().mesh();
//
//             fileName fvPath(mesh.time().path()/"VTK");
//             mkDir(fvPath);
//
//             indirectPrimitivePatch::writeVTK
//             (
//                 fvPath/fileName("uncoveredGgiFaces" + name()),
//                 IndirectList<face>
//                 (
//                     localFaces(),
//                     patchToPatch().uncoveredMasterFaces()
//                 ),
//                 localPoints()
//             );
//         }
//
//         if
//         (
//             !shadow().empty()
//          && patchToPatch().uncoveredSlaveFaces().size() > 0
//         )
//         {
//             // Write uncovered master faces
//             Info<< "Writing uncovered shadow faces for patch "
//                 << shadowName() << " as VTK." << endl;
//
//             const polyMesh& mesh = boundaryMesh().mesh();
//
//             fileName fvPath(mesh.time().path()/"VTK");
//             mkDir(fvPath);
//             Pout<< "Patch " << name()
//                 << " shadow().localFaces(): " << shadow().localFaces().size()
//                 << " patchToPatch().uncoveredSlaveFaces().size(): "
//                 << patchToPatch().uncoveredSlaveFaces().size()
//                 << " shadow().localPoints(): " << shadow().localPoints().size()
//                 << endl;
//
//             indirectPrimitivePatch::writeVTK
//             (
//                 fvPath/fileName("uncoveredGgiFaces" + shadowName()),
//                 IndirectList<face>
//                 (
//                     shadow().localFaces(),
//                     patchToPatch().uncoveredSlaveFaces()
//                 ),
//                 shadow().localPoints()
//             );
//         }

        // Check for bridge overlap
        if (bridgeOverlap())
        {
            InfoInFunction
                << "ggi patch " << name() << " with shadow "
                << shadowName() << " has "
                << patchToPatch().uncoveredMasterFaces().size()
                << " uncovered master faces and "
                << patchToPatch().uncoveredSlaveFaces().size()
                << " uncovered slave faces.  Bridging is switched on. "
                << endl;
        }
    }
}


void Foam::ggiPolyPatch::initOrder(PstreamBuffers&, const primitivePatch&) const
{}


bool Foam::ggiPolyPatch::order
(
    PstreamBuffers& buf,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size(), -1);
    rotation.setSize(pp.size(), 0);

    // Nothing changes
    return false;
}


void Foam::ggiPolyPatch::syncOrder() const
{}


void Foam::ggiPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("shadowPatch") << shadowName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("zone") << zoneName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("bridgeOverlap") << bridgeOverlap_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
