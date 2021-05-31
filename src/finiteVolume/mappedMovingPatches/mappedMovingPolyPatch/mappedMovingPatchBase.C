/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "mappedMovingPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearchMeshObject.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "treeDataFace.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "Time.H"
#include "mapDistribute.H"
#include "SubField.H"
#include "triPointRef.H"
#include "syncTools.H"
#include "treeDataCell.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedMovingPatchBase, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::mappedMovingPatchBase::facePoints
(
    const polyPatch& pp
) const
{
    const polyMesh& mesh = pp.boundaryMesh().mesh();

    // Force construction of min-tet decomp
    (void)mesh.tetBasePtIs();

    // Initialise to face-centre
    tmp<pointField> tfacePoints(new pointField(patch_.size()));
    pointField& facePoints = tfacePoints.ref();

    forAll(pp, facei)
    {
        facePoints[facei] = facePoint
        (
            mesh,
            pp.start()+facei,
            polyMesh::FACE_DIAG_TRIS
        ).rawPoint();
    }

    return tfacePoints;
}


void Foam::mappedMovingPatchBase::collectSamples
(
    const pointField& facePoints,
    pointField& samples,
    labelList& patchFaceProcs,
    labelList& patchFaces,
    pointField& patchFc
) const
{
    // Collect all sample points and the faces they come from.
    {
        List<pointField> globalFc(Pstream::nProcs());
        globalFc[Pstream::myProcNo()] = facePoints;
        Pstream::gatherList(globalFc);
        Pstream::scatterList(globalFc);
        // Rework into straight list
        patchFc = ListListOps::combine<pointField>
        (
            globalFc,
            accessOp<pointField>()
        );
    }

    {
        List<pointField> globalSamples(Pstream::nProcs());
        globalSamples[Pstream::myProcNo()] = facePoints;
        Pstream::gatherList(globalSamples);
        Pstream::scatterList(globalSamples);
        // Rework into straight list
        samples = ListListOps::combine<pointField>
        (
            globalSamples,
            accessOp<pointField>()
        );
    }

    {
        labelListList globalFaces(Pstream::nProcs());
        globalFaces[Pstream::myProcNo()] = identity(patch_.size());
        // Distribute to all processors
        Pstream::gatherList(globalFaces);
        Pstream::scatterList(globalFaces);

        patchFaces = ListListOps::combine<labelList>
        (
            globalFaces,
            accessOp<labelList>()
        );
    }

    {
        labelList nPerProc(Pstream::nProcs());
        nPerProc[Pstream::myProcNo()] = patch_.size();
        Pstream::gatherList(nPerProc);
        Pstream::scatterList(nPerProc);

        patchFaceProcs.setSize(patchFaces.size());

        label sampleI = 0;
        forAll(nPerProc, proci)
        {
            for (label i = 0; i < nPerProc[proci]; i++)
            {
                patchFaceProcs[sampleI++] = proci;
            }
        }
    }
}


// Find the processor/cell containing the samples. Does not account
// for samples being found in two processors.
void Foam::mappedMovingPatchBase::findSamples
(
    const pointField& samples,
    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    // All the info for nearest. Construct to miss
    List<nearInfo> nearest(samples.size());

    const polyPatch& pp = samplePolyPatch();

    if (pp.empty())
    {
        forAll(samples, sampleI)
        {
            nearest[sampleI].second().first() = Foam::sqr(great);
            nearest[sampleI].second().second() = Pstream::myProcNo();
        }
    }
    else
    {
        // Displace face centres by displacement field
        pointField points(pp.faceCentres());
        if (DPtr_)
        {
            points += DPtr_->boundaryField()[pp.index()];
        }

        // Create bound box from
        treeBoundBox patchBb
        (
            treeBoundBox(points).extend(1e-4)
        );

        // If the span is small, use the bounds of the undeformed mesh
        // Happens when only a single face is on the patch
        if (mag(patchBb.span()) < small)
        {
            vector span
            (
                treeBoundBox
                (
                    pp.points(), pp.meshPoints()
                ).extend(1e-4).span()
            );
            patchBb.min() -= span/2.0;
            patchBb.max() += span/2.0;
        }

        indexedOctree<treeDataPoint> boundaryTree
        (
            treeDataPoint(points),
            patchBb,        // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        forAll(samples, sampleI)
        {
            const point& sample = samples[sampleI];

            pointIndexHit& nearInfo = nearest[sampleI].first();
            nearInfo = boundaryTree.findNearest
            (
                sample,
                magSqr(patchBb.span())
            );

            if (!nearInfo.hit())
            {
                nearest[sampleI].second().first() = Foam::sqr(great);
                nearest[sampleI].second().second() =
                    Pstream::myProcNo();
            }
            else
            {
                point fc(points[nearInfo.index()]);
                nearInfo.setPoint(fc);
                nearest[sampleI].second().first() = magSqr(fc-sample);
                nearest[sampleI].second().second() =
                    Pstream::myProcNo();
            }
        }
    }

    // Find nearest. Combine on master.
    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);


    if (debug)
    {
        InfoInFunction
            << "mesh " << sampleRegion() << " : " << endl;

        forAll(nearest, sampleI)
        {
            label proci = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< samples[sampleI]
                << " found on processor:" << proci
                << " in local cell/face/point:" << localI
                << " with location:" << nearest[sampleI].first().rawPoint()
                << endl;
        }
    }

    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());
    sampleLocations.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        if (!nearest[sampleI].first().hit())
        {
            sampleProcs[sampleI] = -1;
            sampleIndices[sampleI] = -1;
            sampleLocations[sampleI] = vector::max;
        }
        else
        {
            sampleProcs[sampleI] = nearest[sampleI].second().second();
            sampleIndices[sampleI] = nearest[sampleI].first().index();
            sampleLocations[sampleI] = nearest[sampleI].first().hitPoint();
        }
    }
}


void Foam::mappedMovingPatchBase::calcMapping() const
{
    if (mapPtr_.valid())
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    // Get points on face (since cannot use face-centres - might be off
    // face-diagonal decomposed tets.
    tmp<pointField> patchPoints(facePoints(patch_));


    // Get global list of all samples and the processor and face they come from.
    pointField samples;
    labelList patchFaceProcs;
    labelList patchFaces;
    pointField patchFc;
    collectSamples
    (
        patchPoints,
        samples,
        patchFaceProcs,
        patchFaces,
        patchFc
    );

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples(samples, sampleProcs, sampleIndices, sampleLocations);

    bool mapSucceeded = true;

    forAll(samples, i)
    {
        if (sampleProcs[i] == -1)
        {
            mapSucceeded = false;
            break;
        }
    }

    if (!mapSucceeded)
    {
        FatalErrorInFunction
            << "Mapping failed for " << nl
            << "    patch: " << patch_.name() << nl
            << "    sampleRegion: " << sampleRegion() << nl
            << "    samplePatch: " << samplePatch() << nl
            << exit(FatalError);
    }

    if (debug && Pstream::master())
    {
        OFstream str
        (
            patch_.boundaryMesh().mesh().time().path()
          / patch_.name()
          + "_mappedMoving.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres/points to file " << str.name()
            << endl;

        label vertI = 0;

        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }

    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchFaceProcs));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.

    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();

    forAll(subMap, proci)
    {
        subMap[proci] = UIndirectList<label>
        (
            sampleIndices,
            subMap[proci]
        );
        constructMap[proci] = UIndirectList<label>
        (
            patchFaces,
            constructMap[proci]
        );

        // if (debug)
        //{
        //    Pout<< "To proc:" << proci << " sending values of cells/faces:"
        //        << subMap[proci] << endl;
        //    Pout<< "From proc:" << proci
        //        << " receiving values of patch faces:"
        //        << constructMap[proci] << endl;
        //}
    }

    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        PackedBoolList used(patch_.size());
        forAll(constructMap, proci)
        {
            const labelList& map = constructMap[proci];

            forAll(map, i)
            {
                label facei = map[i];

                if (used.get(facei) == 0)
                {
                    used.set(facei, 1);
                }
                else
                {
                    FatalErrorInFunction
                        << "On patch " << patch_.name()
                        << " patchface " << facei
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
            }
        }
        forAll(used, facei)
        {
            if (used.get(facei) == 0)
            {
                FatalErrorInFunction
                    << "On patch " << patch_.name()
                    << " patchface " << facei
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp
)
:
    patch_(pp),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    samplePatch_(""),
    coupleGroup_(),
    DPtr_(nullptr),
    mapPtr_(nullptr)
{}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const word& samplePatch
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    samplePatch_(samplePatch),
    coupleGroup_(),
    DPtr_(nullptr),
    mapPtr_(nullptr)
{}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleRegion_(dict.lookupOrDefault<word>("sampleRegion", "")),
    samplePatch_(dict.lookup<word>("samplePatch", "")),
    coupleGroup_(dict),
    DPtr_(nullptr),
    mapPtr_(nullptr)
{
    if (!coupleGroup_.valid())
    {
        if (sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
        }
    }
}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const mappedMovingPatchBase& mpb
)
:
    patch_(pp),
    sampleRegion_(mpb.sampleRegion_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    DPtr_(nullptr),
    mapPtr_(nullptr)
{}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const mappedMovingPatchBase& mpb,
    const labelUList& mapAddressing
)
:
    patch_(pp),
    sampleRegion_(mpb.sampleRegion_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    DPtr_(nullptr),
    mapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingPatchBase::~mappedMovingPatchBase()
{
    clearOut();
}


void Foam::mappedMovingPatchBase::clearOut()
{
    mapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedMovingPatchBase::sampleMesh() const
{
    return patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
    (
        sampleRegion()
    );
}


const Foam::polyPatch& Foam::mappedMovingPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch()
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


Foam::tmp<Foam::pointField> Foam::mappedMovingPatchBase::samplePoints() const
{
    return facePoints(patch_);
}


Foam::pointIndexHit Foam::mappedMovingPatchBase::facePoint
(
    const polyMesh& mesh,
    const label facei,
    const polyMesh::cellDecomposition decompMode
)
{
    const point& fc = mesh.faceCentres()[facei];

    switch (decompMode)
    {
        case polyMesh::FACE_PLANES:
        case polyMesh::FACE_CENTRE_TRIS:
        {
            // For both decompositions the face centre is guaranteed to be
            // on the face
            return pointIndexHit(true, fc, facei);
        }
        break;

        case polyMesh::FACE_DIAG_TRIS:
        case polyMesh::CELL_TETS:
        {
            // Find the intersection of a ray from face centre to cell centre
            // Find intersection of (face-centre-decomposition) centre to
            // cell-centre with face-diagonal-decomposition triangles.

            const pointField& p = mesh.points();
            const face& f = mesh.faces()[facei];

            if (f.size() <= 3)
            {
                // Return centre of triangle.
                return pointIndexHit(true, fc, 0);
            }

            label celli = mesh.faceOwner()[facei];
            const point& cc = mesh.cellCentres()[celli];
            vector d = fc-cc;

            const label fp0 = mesh.tetBasePtIs()[facei];
            const point& basePoint = p[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = p[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = p[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);
                pointHit hitInfo = tri.intersection
                (
                    cc,
                    d,
                    intersection::algorithm::halfRay
                );

                if (hitInfo.hit() && hitInfo.distance() > 0)
                {
                    return pointIndexHit(true, hitInfo.hitPoint(), i-2);
                }

                fp = nextFp;
            }

            // Fall-back
            return pointIndexHit(false, fc, -1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "problem" << abort(FatalError);
            return pointIndexHit();
        }
    }
}


void Foam::mappedMovingPatchBase::write(Ostream& os) const
{
    if (!sampleRegion_.empty())
    {
        writeEntry(os, "sampleRegion", sampleRegion_);
    }
    if (!samplePatch_.empty())
    {
        writeEntry(os, "samplePatch", samplePatch_);
    }

    coupleGroup_.write(os);
}


// ************************************************************************* //
