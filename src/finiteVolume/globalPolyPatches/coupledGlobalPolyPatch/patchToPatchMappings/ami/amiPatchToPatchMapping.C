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

#include "amiPatchToPatchMapping.H"
#include "addToRunTimeSelectionTable.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(amiPatchToPatchMapping, 0);
addToRunTimeSelectionTable
(
    patchToPatchMapping, amiPatchToPatchMapping, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void amiPatchToPatchMapping::makeInterpolator() const
{
    if (interpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer is already set!"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "Create AMI zone-to-zone interpolator" << endl;
    }

    interpolatorPtr_.set
    (
        new amiZoneInterpolation
        (
            zoneA(),
            zoneB(),
            faceAreaIntersect::tmMesh, // triMode
            true,   // requireMatch
            amiZoneInterpolation::imFaceAreaWeight, // interpolationMethodNames
            -1,     // lowWeightCorrection
            false,  // reverseTarget
            true    // use globalPolyPatch
        )
    );

    checkZoneAToZoneBError();
    //checkZoneBToZoneAError();
}


const amiZoneInterpolation&
amiPatchToPatchMapping::interpolator() const
{
    if (interpolatorPtr_.empty())
    {
        makeInterpolator();
    }

    return interpolatorPtr_();
}


void amiPatchToPatchMapping::checkZoneAToZoneBError() const
{
    // Reference to patch face centres
    const vectorField& patchAFaceCentres = patchA().faceCentres();

    // Construct global zone field
    const vectorField zoneAFaceCentres
    (
        globalPatchA().patchFaceToGlobal(patchAFaceCentres)
    );

    // Interpolate global zone field from A to B
    const vectorField zoneBFaceCentres
    (
        interpolator().interpolateToTarget(zoneAFaceCentres)
    );

    // Extract local patch field
    const vectorField patchBFaceCentres
    (
        globalPatchB().globalFaceToPatch(zoneBFaceCentres)
    );

    // Print maximum error
    if (debug)
    {
        Info<< "interface-to-interface face error: "
            << gMax(mag(patchBFaceCentres - patchB().faceCentres()))
            << endl;
    }

}


void amiPatchToPatchMapping::calcZoneAPointAddressing() const
{
    if (zoneAPointAddressingPtr_)
    {
        FatalErrorIn(type())
            << "zoneA points addressing already exists"
            << abort(FatalError);
    }

    zoneAPointAddressingPtr_ =
        new List<labelPair>
        (
            zoneA().nPoints(), labelPair(-1,-1)
        );
    List<labelPair>& zoneAPointAddr = *zoneAPointAddressingPtr_;

    const labelListList& zoneAFaceAddr = interpolator().srcAddress();
    const labelListList& pointFaces = zoneA().pointFaces();
    const faceList& zoneBFaces = zoneB().localFaces();
    const pointField& zoneBPoints = zoneB().localPoints();
    const pointField& zoneAPoints = zoneA().localPoints();

    forAll(zoneAPointAddr, pointI)
    {
        const point& P = zoneAPoints[pointI];
        labelHashSet possibleZoneBFacesSet;
        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            const label curFace = curPointFaces[faceI];
            const labelList& curZoneBFaces = zoneAFaceAddr[curFace];

            forAll(curZoneBFaces, fI)
            {
                if (!possibleZoneBFacesSet.found(curZoneBFaces[fI]))
                {
                    possibleZoneBFacesSet.insert(curZoneBFaces[fI]);
                }
            }
        }

        const labelList possibleZoneBFaces = possibleZoneBFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);

        forAll(possibleZoneBFaces, faceI)
        {
            const label curZoneBFace = possibleZoneBFaces[faceI];
            const face& f = zoneBFaces[curZoneBFace];
            const point ctr = Foam::average(f.points(zoneBPoints));
            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = zoneBPoints[f.nextLabel(pI)];

                const triPointRef t
                (
                    zoneBPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                const scalar A = mag(n);
                n /= A;

                // Intersection point
                const point I = P + n*(n & (t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                const scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curZoneBFace;
                    faceTriangle.second() = pI;

                    //distance = ((P - I) & n);
                }
            }
        }

        zoneAPointAddr[pointI] = faceTriangle;
    }


    // Check orientation

    const pointField& zoneAPointNormals = zoneA().pointNormals();

    const vectorField& zoneBFaceNormals = zoneB().faceNormals();

    scalarField orientation(zoneAPointAddr.size(), 0);

    label nIncorrectPoints = 0;

    forAll(zoneAPointAddr, pointI)
    {
        orientation[pointI] =
            (
                zoneAPointNormals[pointI]
              & zoneBFaceNormals[zoneAPointAddr[pointI].first()]
            );

        if (orientation[pointI] > -SMALL)
        {
            nIncorrectPoints++;
        }
    }

    if (debug)
    {
        Info<< "zoneA point orientation (< 0), max: "
            << max(orientation)
            << ", min: " << min(orientation) << ", nIncorrectPoints: "
            << nIncorrectPoints << "/" << zoneAPointAddr.size() << endl;
    }
}


void amiPatchToPatchMapping::calcZoneAPointWeights() const
{
    if (zoneAPointWeightsPtr_)
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    zoneAPointWeightsPtr_ =
        new FieldField<Field, scalar>(zoneA().nPoints());
    FieldField<Field, scalar>& zoneAPointWeights = *zoneAPointWeightsPtr_;

    const faceList& zoneBFaces = zoneB().localFaces();
    const pointField& zoneBPoints = zoneB().localPoints();
    const pointField& zoneAPoints = zoneA().localPoints();

    const List<labelPair>& addr = this->zoneAPointAddr();

    forAll(zoneAPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = zoneAPoints[pointI];
            const face& hitFace = zoneBFaces[addr[pointI].first()];
            const point ctr = Foam::average(hitFace.points(zoneBPoints));
            const label pI = addr[pointI].second();

            const triPointRef t
            (
                zoneBPoints[hitFace[pI]],
                zoneBPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            const point I = P + n*(n&(t.a() - P));

            zoneAPointWeights.set(pointI, new scalarField(3));

            zoneAPointWeights[pointI][0] = t.pointToBarycentric(I).a();
            zoneAPointWeights[pointI][1] = t.pointToBarycentric(I).b();
            zoneAPointWeights[pointI][2] = t.pointToBarycentric(I).c();
        }
        else
        {
            zoneAPointWeights.set(pointI, new scalarField(0));
        }
    }
}


void amiPatchToPatchMapping::calcZoneBPointAddressing() const
{
    if (zoneBPointAddressingPtr_)
    {
        FatalErrorIn(type())
            << "zoneB points addressing already exists"
            << abort(FatalError);
    }

    zoneBPointAddressingPtr_ =
        new List<labelPair>
        (
            zoneB().nPoints(), labelPair(-1,-1)
        );
    List<labelPair>& zoneBPointAddr = *zoneBPointAddressingPtr_;

    const labelListList& zoneBFaceAddr = interpolator().tgtAddress();
    const labelListList& pointFaces = zoneB().pointFaces();
    const faceList& zoneAFaces = zoneA().localFaces();
    const pointField& zoneAPoints = zoneA().localPoints();
    const pointField& zoneBPoints = zoneB().localPoints();

    forAll(zoneBPointAddr, pointI)
    {
        const point& P = zoneBPoints[pointI];
        labelHashSet possibleZoneAFacesSet;
        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            const label curFace = curPointFaces[faceI];
            const labelList& curZoneAFaces = zoneBFaceAddr[curFace];

            forAll(curZoneAFaces, fI)
            {
                if (!possibleZoneAFacesSet.found(curZoneAFaces[fI]))
                {
                    possibleZoneAFacesSet.insert(curZoneAFaces[fI]);
                }
            }
        }

        const labelList possibleZoneAFaces = possibleZoneAFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);

        forAll(possibleZoneAFaces, faceI)
        {
            const label curZoneAFace = possibleZoneAFaces[faceI];
            const face& f = zoneAFaces[curZoneAFace];
            const point ctr = Foam::average(f.points(zoneAPoints));
            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = zoneAPoints[f.nextLabel(pI)];

                const triPointRef t
                (
                    zoneAPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                const scalar A = mag(n);
                n /= A;

                // Intersection point
                const point I = P + n*(n & (t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                const scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curZoneAFace;
                    faceTriangle.second() = pI;

                    //distance = ((P - I) & n);
                }
            }
        }

        zoneBPointAddr[pointI] = faceTriangle;
    }

    // Check orientation
    const pointField& zoneBPointNormals = zoneB().pointNormals();
    const vectorField& zoneAFaceNormals = zoneA().faceNormals();

    scalarField orientation(zoneBPointAddr.size(), 0);

    label nIncorrectPoints = 0;

    forAll(zoneBPointAddr, pointI)
    {
        orientation[pointI] =
            (
                zoneBPointNormals[pointI]
              & zoneAFaceNormals[zoneBPointAddr[pointI].first()]
            );

        if (orientation[pointI] > -SMALL)
        {
            nIncorrectPoints++;
        }
    }

    if (debug)
    {
        Info<< "zoneB point orientation (< 0), max: "
            << max(orientation)
            << ", min: " << min(orientation) << ", nIncorrectPoints: "
            << nIncorrectPoints << "/" << zoneBPointAddr.size() << endl;
    }
}


void amiPatchToPatchMapping::calcZoneBPointWeights() const
{
    if (zoneBPointWeightsPtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }
    zoneBPointWeightsPtr_ =
        new FieldField<Field, scalar>(zoneB().nPoints());
    FieldField<Field, scalar>& zoneBPointWeights = *zoneBPointWeightsPtr_;

    const faceList& zoneAFaces = zoneA().localFaces();
    const pointField& zoneAPoints = zoneA().localPoints();
    const pointField& zoneBPoints = zoneB().localPoints();
    const List<labelPair>& addr = this->zoneBPointAddr();

    forAll(zoneBPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = zoneBPoints[pointI];
            const face& hitFace = zoneAFaces[addr[pointI].first()];
            const point ctr = Foam::average(hitFace.points(zoneAPoints));
            const label pI = addr[pointI].second();

            const triPointRef t
            (
                zoneAPoints[hitFace[pI]],
                zoneAPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            const point I = P + n*(n&(t.a() - P));

            zoneBPointWeights.set(pointI, new scalarField(3));

            zoneBPointWeights[pointI][0] = t.pointToBarycentric(I).a();
            zoneBPointWeights[pointI][1] = t.pointToBarycentric(I).b();
            zoneBPointWeights[pointI][2] = t.pointToBarycentric(I).c();
        }
        else
        {
            zoneBPointWeights.set(pointI, new scalarField(0));
        }
    }
}

const List<labelPair>&
amiPatchToPatchMapping::zoneAPointAddr() const
{
    if (!zoneAPointAddressingPtr_)
    {
        calcZoneAPointAddressing();
    }

    return *zoneAPointAddressingPtr_;
}


const FieldField<Field, scalar>&
amiPatchToPatchMapping::zoneAPointWeights() const
{
    if (!zoneAPointWeightsPtr_)
    {
        calcZoneAPointWeights();
    }

    return *zoneAPointWeightsPtr_;
}


const List<labelPair>&
amiPatchToPatchMapping::zoneBPointAddr() const
{
    if (!zoneBPointAddressingPtr_)
    {
        calcZoneBPointAddressing();
    }

    return *zoneBPointAddressingPtr_;
}


const FieldField<Field, scalar>&
amiPatchToPatchMapping::zoneBPointWeights() const
{
    if (!zoneBPointWeightsPtr_)
    {
        calcZoneBPointWeights();
    }

    return *zoneBPointWeightsPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

amiPatchToPatchMapping::amiPatchToPatchMapping
(
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalPatchA,
    const globalPolyPatch& globalPatchB
)
:
    patchToPatchMapping
    (
        typeName_(), dict, patchA, patchB, globalPatchA, globalPatchB
    ),
    interpolatorPtr_(nullptr),
    zoneAPointAddressingPtr_(nullptr),
    zoneAPointWeightsPtr_(nullptr),
    zoneBPointAddressingPtr_(nullptr),
    zoneBPointWeightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void amiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFaces<scalar>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFaces<vector>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferFaces<symmTensor>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferFaces<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferFaces<tensor>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPoints<scalar>(fromZone, toZone, fromField, toField);
}

void amiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferPoints<vector>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferPoints<symmTensor>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferPoints<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void amiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferPoints<tensor>(fromZone, toZone, fromField, toField);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings
} // End namespace Foam

// ************************************************************************* //
