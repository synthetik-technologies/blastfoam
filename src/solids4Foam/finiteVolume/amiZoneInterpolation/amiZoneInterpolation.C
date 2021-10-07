/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#ifdef OPENFOAMESIORFOUNDATION

#include "amiZoneInterpolation.H"
#ifdef OPENFOAMESI
    #include "AMIMethod.H"
    #include "directAMI.H"
    #include "mapNearestAMI.H"
    #include "faceAreaWeightAMI.H"
    #include "partialFaceAreaWeightAMI.H"
#else
    #include "newAMIMethod.H"
    #include "newDirectAMI.H"
    #include "newMapNearestAMI.H"
    #include "newFaceAreaWeightAMI.H"
    #include "newPartialFaceAreaWeightAMI.H"
    #include "newSweptFaceAreaWeightAMI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(amiZoneInterpolation, 0);

#ifdef OPENFOAMESI
    makeAMIMethod(amiZoneInterpolation);
    makeAMIMethodType(amiZoneInterpolation, directAMI);
    makeAMIMethodType(amiZoneInterpolation, mapNearestAMI);
    makeAMIMethodType(amiZoneInterpolation, faceAreaWeightAMI);
    makeAMIMethodType(amiZoneInterpolation, partialFaceAreaWeightAMI);
#else
    makeAMIMethod(amiZoneInterpolation);
    makeAMIMethodType(amiZoneInterpolation, newDirectAMI);
    makeAMIMethodType(amiZoneInterpolation, newMapNearestAMI);
    makeAMIMethodType(amiZoneInterpolation, newFaceAreaWeightAMI);
    makeAMIMethodType(amiZoneInterpolation, newPartialFaceAreaWeightAMI);
    makeAMIMethodType(amiZoneInterpolation, newSweptFaceAreaWeightAMI);
#endif
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::amiZoneInterpolation::calcSourcePointAddressing() const
{
    // Find source points addressing
    if (sourcePointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void amiZoneInterpolation::"
            "calcSourcePointAddressing() const"
        )
            << "Source points addressing already exists"
                << abort(FatalError);
    }

    sourcePointAddressingPtr_ =
        new List<labelPair>
        (
            sourcePatch_.nPoints(),
            labelPair(-1,-1)
        );
    List<labelPair>& sourcePointAddr = *sourcePointAddressingPtr_;

    sourcePointDistancePtr_ =
        new scalarField
        (
            sourcePatch_.nPoints(),
            GREAT
        );
    scalarField& sourcePointDist = *sourcePointDistancePtr_;

    const labelListList& sourceFaceAddr = this->srcAddress();

    const labelListList& pointFaces = sourcePatch_.pointFaces();

    const faceList& targetFaces = targetPatch_.localFaces();
    const pointField& targetPoints = targetPatch_.localPoints();

    const pointField& sourcePoints = sourcePatch_.localPoints();

    forAll(sourcePointAddr, pointI)
    {
        const point& P = sourcePoints[pointI];

        labelHashSet possibleTargetFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            const labelList& curTargetFaces = sourceFaceAddr[curFace];

            forAll(curTargetFaces, fI)
            {
                if (!possibleTargetFacesSet.found(curTargetFaces[fI]))
                {
                    possibleTargetFacesSet.insert(curTargetFaces[fI]);
                }
            }
        }

        labelList possibleTargetFaces = possibleTargetFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleTargetFaces, faceI)
        {
            label curTargetFace = possibleTargetFaces[faceI];

            const face& f = targetFaces[curTargetFace];

            point ctr = Foam::average(f.points(targetPoints));

            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = targetPoints[f.nextLabel(pI)];

                triPointRef t
                (
                    targetPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                scalar A = mag(n);
                n /= A;

                // Intersection point
                point I = P + n*(n&(t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curTargetFace;
                    faceTriangle.second() = pI;

                    distance = ((P - I)&n);
                }
            }
        }

        sourcePointAddr[pointI] = faceTriangle;
        sourcePointDist[pointI] = distance;
    }
}


void Foam::amiZoneInterpolation::calcSourcePointWeights() const
{
    // Find source point weights
    if (sourcePointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void amiZoneInterpolation::"
            "calcSourcePointAddressing() const"
        )
            << "Source point weights already exist"
                << abort(FatalError);
    }

    // sourcePointWeightsPtr_ =
    //     new FieldField<Field, scalar>(sourcePatch_.nPoints());
    // FieldField<Field, scalar>& sourcePointWeights = *sourcePointWeightsPtr_;
    sourcePointWeightsPtr_ =
        new List< List<scalar> >(sourcePatch_.nPoints());
    List< List<scalar> >& sourcePointWeights = *sourcePointWeightsPtr_;

    const faceList& targetFaces = targetPatch_.localFaces();
    const pointField& targetPoints = targetPatch_.localPoints();

    const pointField& sourcePoints = sourcePatch_.localPoints();

    const List<labelPair>& addr = this->sourcePointAddr();

    forAll(sourcePointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = sourcePoints[pointI];

            const face& hitFace =
                targetFaces[addr[pointI].first()];

            point ctr = Foam::average(hitFace.points(targetPoints));

            label pI = addr[pointI].second();

            triPointRef t
            (
                targetPoints[hitFace[pI]],
                targetPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            point I = P + n*(n&(t.a() - P));

            sourcePointWeights[pointI] = List<scalar>(3);

            // sourcePointWeights[pointI][0] = t.Ni(0, I);
            // sourcePointWeights[pointI][1] = t.Ni(1, I);
            // sourcePointWeights[pointI][2] = t.Ni(2, I);

            sourcePointWeights[pointI][0] = t.pointToBarycentric(I).a();
            sourcePointWeights[pointI][1] = t.pointToBarycentric(I).b();
            sourcePointWeights[pointI][2] = t.pointToBarycentric(I).c();
        }
        else
        {
            sourcePointWeights[pointI] = List<scalar>(0);
        }
    }
}


void Foam::amiZoneInterpolation::calcTargetPointAddressing() const
{
    // Find source points addressing
    if (targetPointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void amiZoneInterpolation::"
            "calcTargetPointAddressing() const"
        )
            << "Target points addressing already exists"
                << abort(FatalError);
    }

    targetPointAddressingPtr_ =
        new List<labelPair>
        (
            targetPatch_.nPoints(),
            labelPair(-1,-1)
        );
    List<labelPair>& targetPointAddr = *targetPointAddressingPtr_;

    targetPointDistancePtr_ =
        new scalarField
        (
            targetPatch_.nPoints(),
            GREAT
        );
    scalarField& targetPointDist = *targetPointDistancePtr_;

    const labelListList& targetFaceAddr = this->tgtAddress();

    const labelListList& pointFaces = targetPatch_.pointFaces();

    const faceList& sourceFaces = sourcePatch_.localFaces();
    const pointField& sourcePoints = sourcePatch_.localPoints();

    const pointField& targetPoints = targetPatch_.localPoints();

    forAll(targetPointAddr, pointI)
    {
        const point& P = targetPoints[pointI];

        labelHashSet possibleSourceFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            const labelList& curSourceFaces = targetFaceAddr[curFace];

            forAll(curSourceFaces, fI)
            {
                if (!possibleSourceFacesSet.found(curSourceFaces[fI]))
                {
                    possibleSourceFacesSet.insert(curSourceFaces[fI]);
                }
            }
        }

        labelList possibleSourceFaces = possibleSourceFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleSourceFaces, faceI)
        {
            label curSourceFace = possibleSourceFaces[faceI];

            const face& f = sourceFaces[curSourceFace];

            point ctr = Foam::average(f.points(sourcePoints));

            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = sourcePoints[f.nextLabel(pI)];

                triPointRef t
                (
                    sourcePoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                scalar A = mag(n);
                n /= A;

                // Intersection point
                point I = P + n*(n&(t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curSourceFace;
                    faceTriangle.second() = pI;

                    distance = ((P - I)&n);
                }
            }
        }

//         Info << "MinEta " << MinEta << endl;

        targetPointAddr[pointI] = faceTriangle;
        targetPointDist[pointI] = distance;
    }
}


void Foam::amiZoneInterpolation::calcTargetPointWeights() const
{
    // Find source point weights
    if (targetPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void amiZoneInterpolation::"
            "calcTargetPointAddressing() const"
        )
            << "Target point weights already exist"
                << abort(FatalError);
    }

    // targetPointWeightsPtr_ =
    //     new FieldField<Field, scalar>(targetPatch_.nPoints());
    // FieldField<Field, scalar>& targetPointWeights = *targetPointWeightsPtr_;
    targetPointWeightsPtr_ =
        new List< List<scalar> >(targetPatch_.nPoints());
    List< List<scalar> >& targetPointWeights = *targetPointWeightsPtr_;

    const faceList& sourceFaces = sourcePatch_.localFaces();
    const pointField& sourcePoints = sourcePatch_.localPoints();

    const pointField& targetPoints = targetPatch_.localPoints();

    const List<labelPair>& addr = this->targetPointAddr();

    forAll(targetPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = targetPoints[pointI];

            const face& hitFace =
                sourceFaces[addr[pointI].first()];

            point ctr = Foam::average(hitFace.points(sourcePoints));

            label pI = addr[pointI].second();

            triPointRef t
            (
                sourcePoints[hitFace[pI]],
                sourcePoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            point I = P + n*(n&(t.a() - P));

            targetPointWeights[pointI] = List<scalar>(3);

            // targetPointWeights[pointI][0] = t.Ni(0, I);
            // targetPointWeights[pointI][1] = t.Ni(1, I);
            // targetPointWeights[pointI][2] = t.Ni(2, I);

            targetPointWeights[pointI][0] = t.pointToBarycentric(I).a();
            targetPointWeights[pointI][1] = t.pointToBarycentric(I).b();
            targetPointWeights[pointI][2] = t.pointToBarycentric(I).c();
        }
        else
        {
            targetPointWeights[pointI] = List<scalar>(0);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::amiZoneInterpolation::amiZoneInterpolation
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const bool useGlobalPolyPatch
)
:
#ifdef OPENFOAMESI
    AMIInterpolation<standAlonePatch, standAlonePatch>
#else
    newAMIInterpolation<standAlonePatch, standAlonePatch>
#endif
    (
        srcPatch,
        tgtPatch,
        triMode,
        true,
        imFaceAreaWeight,
        -1,
        false,
        useGlobalPolyPatch
    ),
    sourcePatch_(srcPatch),
    targetPatch_(tgtPatch),
    sourcePatchInterp_(sourcePatch_),
    targetPatchInterp_(targetPatch_),
    sourcePointAddressingPtr_(nullptr),
    sourcePointWeightsPtr_(nullptr),
    sourcePointDistancePtr_(nullptr),
    targetPointAddressingPtr_(nullptr),
    targetPointWeightsPtr_(nullptr),
    targetPointDistancePtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::amiZoneInterpolation::~amiZoneInterpolation()
{
    deleteDemandDrivenData(sourcePointAddressingPtr_);
    deleteDemandDrivenData(sourcePointWeightsPtr_);
    deleteDemandDrivenData(sourcePointDistancePtr_);
    deleteDemandDrivenData(targetPointAddressingPtr_);
    deleteDemandDrivenData(targetPointWeightsPtr_);
    deleteDemandDrivenData(targetPointDistancePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::labelPair>&
Foam::amiZoneInterpolation::sourcePointAddr() const
{
    if (!sourcePointAddressingPtr_)
    {
        calcSourcePointAddressing();
    }

    return *sourcePointAddressingPtr_;
}


const Foam::List< Foam::List<Foam::scalar> >&
Foam::amiZoneInterpolation::sourcePointWeights() const
{
    if (!sourcePointWeightsPtr_)
    {
        calcSourcePointWeights();
    }

    return *sourcePointWeightsPtr_;
}


const Foam::scalarField&
Foam::amiZoneInterpolation::sourcePointDistanceToIntersection() const
{
    if (!sourcePointDistancePtr_)
    {
        calcSourcePointAddressing();
    }

    return *sourcePointDistancePtr_;
}


const Foam::List<Foam::labelPair>&
Foam::amiZoneInterpolation::targetPointAddr() const
{
    if (!targetPointAddressingPtr_)
    {
        calcTargetPointAddressing();
    }

    return *targetPointAddressingPtr_;
}


const Foam::List< Foam::List<Foam::scalar> >&
Foam::amiZoneInterpolation::targetPointWeights() const
{
    if (!targetPointWeightsPtr_)
    {
        calcTargetPointWeights();
    }

    return *targetPointWeightsPtr_;
}


const Foam::scalarField&
Foam::amiZoneInterpolation::targetPointDistanceToIntersection() const
{
    if (!targetPointDistancePtr_)
    {
        calcTargetPointAddressing();
    }

    return *targetPointDistancePtr_;
}

#endif // end of #ifdef OPENFOAMESIORFOUNDATION

// ************************************************************************* //
