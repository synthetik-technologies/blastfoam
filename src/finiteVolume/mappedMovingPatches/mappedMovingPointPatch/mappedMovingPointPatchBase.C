/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technology: Mapping of point patches
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

#include "mappedMovingPointPatchBase.H"
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
    defineTypeNameAndDebug(mappedMovingPointPatchBase, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mappedMovingPointPatchBase::collectSamples
(
    const pointField& points,
    pointField& samples,
    labelList& patchPointProcs,
    labelList& patchPoints,
    pointField& patchP
) const
{
    // Collect all sample points and the faces they come from.
    {
        List<pointField> globalP(Pstream::nProcs());
        globalP[Pstream::myProcNo()] = points;
        Pstream::gatherList(globalP);
        Pstream::scatterList(globalP);
        // Rework into straight list
        patchP = ListListOps::combine<pointField>
        (
            globalP,
            accessOp<pointField>()
        );
    }

    {
        List<pointField> globalSamples(Pstream::nProcs());
        globalSamples[Pstream::myProcNo()] = points;
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
        labelListList globalPoints(Pstream::nProcs());
        globalPoints[Pstream::myProcNo()] = identity(patchPtr_->size());
        // Distribute to all processors
        Pstream::gatherList(globalPoints);
        Pstream::scatterList(globalPoints);

        patchPoints = ListListOps::combine<labelList>
        (
            globalPoints,
            accessOp<labelList>()
        );
    }

    {
        labelList nPerProc(Pstream::nProcs());
        nPerProc[Pstream::myProcNo()] = patchPtr_->size();
        Pstream::gatherList(nPerProc);
        Pstream::scatterList(nPerProc);

        patchPointProcs.setSize(patchPoints.size());

        label sampleI = 0;
        forAll(nPerProc, proci)
        {
            for (label i = 0; i < nPerProc[proci]; i++)
            {
                patchPointProcs[sampleI++] = proci;
            }
        }
    }
}


// Find the processor/cell containing the samples. Does not account
// for samples being found in two processors.
void Foam::mappedMovingPointPatchBase::findSamples
(
    const pointField& samples,
    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    // All the info for nearest. Construct to miss
    List<mappedMovingPatchBase::nearInfo> nearest(samples.size());

    const polyPatch& pp = samplePolyPatch();

    pointField points(pp.localPoints());
    if (displacementPtr_.valid())
    {
        points +=
            displacementPtr_->boundaryField()[pp.index()].patchInternalField();
    }

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
        // patch (local) points
        treeBoundBox patchBb
        (
            treeBoundBox(points).extend(1e-4)
        );

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
                const point& pt = nearInfo.hitPoint();

                nearest[sampleI].second().first() = magSqr(pt-sample);
                nearest[sampleI].second().second() =
                    Pstream::myProcNo();
            }
        }
    }


    // Find nearest. Combine on master.
    Pstream::listCombineGather
    (
        nearest,
        mappedMovingPatchBase::nearestEqOp()
    );
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


void Foam::mappedMovingPointPatchBase::calcMapping() const
{
    if (!patchPtr_.valid())
    {
        patchPtr_.reset
        (
            &pMesh_.lookupObject<pointMesh>
            (
                "pointMesh"
            ).boundary()[pp_.index()]
        );
    }

    if (mapPtr_.valid())
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    // Get points on face (since cannot use face-centres - might be off
    // face-diagonal decomposed tets.
    tmp<pointField> points(patchPtr_->localPoints());

    // Get global list of all samples and the processor and face they come from.
    pointField samples;
    labelList patchPointProcs;
    labelList patchPoints;
    pointField patchP;
    collectSamples
    (
        points,
        samples,
        patchPointProcs,
        patchPoints,
        patchP
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
            << "    patch: " << patchPtr_->name() << nl
            << "    sampleRegion: " << sampleRegion() << nl
            << "    samplePatch: " << samplePatch() << nl
            << exit(FatalError);
    }

    if (debug && Pstream::master())
    {
        OFstream str
        (
            patchPtr_->boundaryMesh().mesh().time().path()
          / patchPtr_->name()
          + "_mappedMoving.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres/points to file " << str.name()
            << endl;

        label vertI = 0;

        forAll(patchP, i)
        {
            meshTools::writeOBJ(str, patchP[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }

    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchPointProcs));

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
            patchPoints,
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
    mapPtr_().constructSize() = patchPtr_->size();

    if (debug)
    {
        // Check that all elements get a value.
        PackedBoolList used(patchPtr_->size());
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
                        << "On patch " << patchPtr_->name()
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
                    << "On patch " << patchPtr_->name()
                    << " patchface " << facei
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }
}


// Hack to read old (List-based) format. See Field.C. The difference
// is only that in case of falling back to old format it expects a non-uniform
// list instead of a single vector.
Foam::tmp<Foam::pointField> Foam::mappedMovingPointPatchBase::readListOrField
(
    const word& keyword,
    const dictionary& dict,
    const label size
)
{
    tmp<pointField> tfld(new pointField());
    pointField& fld = tfld.ref();

    if (size)
    {
        ITstream& is = dict.lookup(keyword);

        // Read first token
        token firstToken(is);

        if (firstToken.isWord())
        {
            if (firstToken.wordToken() == "uniform")
            {
                fld.setSize(size);
                fld = pTraits<vector>(is);
            }
            else if (firstToken.wordToken() == "nonuniform")
            {
                is >> static_cast<List<vector>&>(fld);
                if (fld.size() != size)
                {
                    FatalIOErrorInFunction
                    (
                        dict
                    )   << "size " << fld.size()
                        << " is not equal to the given value of " << size
                        << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', found "
                    << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else
        {
            if (is.version() == 2.0)
            {
                IOWarningInFunction
                (
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', "
                       "assuming List format for backwards compatibility."
                       "Foam version 2.0." << endl;

                is.putBack(firstToken);
                is >> static_cast<List<vector>&>(fld);
            }
        }
    }
    return tfld;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingPointPatchBase::mappedMovingPointPatchBase
(
    const polyPatch& pp
)
:
    pp_(pp),
    mpp_(dynamic_cast<const mappedMovingPatchBase&>(pp)),
    pMesh_(pp.boundaryMesh().mesh()),
    patchPtr_(nullptr),
    displacementPtr_(nullptr),
    mapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingPointPatchBase::~mappedMovingPointPatchBase()
{
    clearOut();
}


void Foam::mappedMovingPointPatchBase::clearOut()
{
    mapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedMovingPointPatchBase::sampleMesh() const
{
    return pMesh_.time().lookupObject<polyMesh>
    (
        sampleRegion()
    );
}


const Foam::polyPatch& Foam::mappedMovingPointPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch()
            << " in region " << sampleRegion() << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


Foam::tmp<Foam::pointField> Foam::mappedMovingPointPatchBase::samplePoints() const
{
    if (!patchPtr_.valid())
    {
        patchPtr_.reset
        (
            &pMesh_.lookupObject<pointMesh>
            (
                "pointMesh"
            ).boundary()[pp_.index()]
        );
    }

    return patchPtr_->localPoints();
}


void Foam::mappedMovingPointPatchBase::write(Ostream& os) const
{}


// ************************************************************************* //
