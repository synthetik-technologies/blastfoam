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

#include "inverseDistancePatchToPatchMapping.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{
    defineTypeNameAndDebug(inverseDistance, 0);
    addToRunTimeSelectionTable
    (
        patchToPatchMapping,
        inverseDistance,
        dictionary
    );
}
}

bool Foam::patchToPatchMappings::inverseDistance::rayHitsFace
(
    const point& p,
    const vector& r,
    const face& f,
    const pointField& ps
)
{
    using namespace constant::mathematical;

    const tensor T = tensor::I - sqr(r);

    scalar angle = 0;

    forAll(f, i)
    {
        const vector& a = T & (ps[f[i]] - p);
        const vector& b = T & (ps[f[f.fcIndex(i)]] - p);

        const scalar meanMagSqrAB = (magSqr(a) + magSqr(b))/2;
        const scalar geometricMeanMagSqrAB = sqrt(magSqr(a)*magSqr(b));

        // This indicates that we have hit a point to within round off error
        if (geometricMeanMagSqrAB < small*meanMagSqrAB) return true;

        angle -=
            sign(r & (a ^ b))
           *acos(min(max(-1, (a & b)/geometricMeanMagSqrAB), +1));
    }

    return pi < angle && angle < 3*pi;
}


void Foam::patchToPatchMappings::inverseDistance::generateFaceWeights
(
    const primitivePatch& patch,
    const primitivePatch& otherPatch,
    const bool reverse,
    List<DynamicList<label>>& otherFaces,
    List<DynamicList<scalar>>& weights
)
{
    forAll(otherFaces, facei)
    {
        if (otherFaces[facei].empty()) continue;

        label otherFacei = -1;

        // Find the other face that "contains" this face's centre
        forAll(otherFaces[facei], i)
        {
            if
            (
                rayHitsFace
                (
                    patch.faceCentres()[facei],
                    (reverse ? -1 : +1)*patch.faceNormals()[facei],
                    otherPatch[otherFaces[facei][i]],
                    otherPatch.points()
                )
            )
            {
                otherFacei = otherFaces[facei][i];
                break;
            }
        }

        const point& c = patch.faceCentres()[facei];

        // If the above failed, find the closest
        if (otherFacei == -1)
        {
            scalar minDistSqr = vGreat;

            forAll(otherFaces[facei], i)
            {
                const point& otherC =
                    otherPatch.faceCentres()[otherFaces[facei][i]];
                const scalar distSqr = magSqr(c - otherC);
                if (distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    otherFacei = otherFaces[facei][i];
                }
            }
        }

        // Remove all faces
        otherFaces[facei].clear();

        // Add the found face and all its neighbours
        otherFaces[facei].append(otherFacei);
        weights[facei].append
        (
            1.0/(mag(c - otherPatch.faceCentres()[otherFacei]) + rootVSmall)
        );

        forAll(otherPatch.faceFaces()[otherFacei], i)
        {
            const label otherFacej = otherPatch.faceFaces()[otherFacei][i];

            otherFaces[facei].append(otherFacej);
            weights[facei].append
            (
                1.0
               /(mag(c - otherPatch.faceCentres()[otherFacej]) + rootVSmall)
            );
        }
    }
}


void Foam::patchToPatchMappings::inverseDistance::generatePointWeights
(
    const primitivePatch& patch,
    const primitivePatch& origOtherPatch,
    const bool merge,
    List<DynamicList<label>>& otherPoints,
    List<DynamicList<scalar>>& weights
)
{
    tmpNrc<primitivePatch> totherPatch(origOtherPatch);
    autoPtr<pointField> otherPointsPtr;
    autoPtr<faceList> otherFacesPtr;
    labelList mergePointMap(identity(origOtherPatch.nPoints()));

    if (merge)
    {
        if
        (
            mergePoints
            (
                origOtherPatch.points(),
                1e-6,
                false,
                mergePointMap
            )
        )
        {
            otherPointsPtr.reset
            (
                new pointField(origOtherPatch.points(), mergePointMap)
            );
            otherFacesPtr.reset(new faceList(origOtherPatch));
            faceList& otherFaces = otherFacesPtr();
            forAll(otherFaces, facei)
            {
                face& f = otherFaces[facei];
                forAll(f, i)
                {
                    f[i] = mergePointMap[f[i]];
                }
            }
            totherPatch =
                new primitivePatch
                (
                    SubList<face>(otherFacesPtr(), otherFacesPtr->size()),
                    otherPointsPtr()
                );
        }
    }
    const primitivePatch& otherPatch = totherPatch();

    Map<label> rMergePointMap;
    forAll(mergePointMap, i)
    {
        rMergePointMap.insert(mergePointMap[i], i);
    }

    const edgeList& otherEdges = otherPatch.edges();
    const labelListList& otherPointEdges = otherPatch.pointEdges();

    forAll(otherPoints, pointi)
    {
        if (otherPoints[pointi].empty()) continue;

        label otherPointi = -1;

        const point& p = patch.points()[pointi];

        scalar minDistSqr = vGreat;

        forAll(otherPoints[pointi], i)
        {
            const point& otherP =
                otherPatch.points()[otherPoints[pointi][i]];
            const scalar distSqr = magSqr(p - otherP);
            if (distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                otherPointi = mergePointMap[otherPoints[pointi][i]];
            }
        }

        // Remove all faces
        otherPoints[pointi].clear();

        // Add the found face and all its neighbours
        const point& otherP = otherPatch.points()[otherPointi];
        otherPoints[pointi].append(rMergePointMap[otherPointi]);
        weights[pointi].append(1/(mag(p - otherP) + rootVSmall));

        const labelList& pEdges = otherPointEdges[otherPointi];
        forAll(pEdges, i)
        {
            const label otherPointj =
                otherEdges[pEdges[i]].otherVertex(otherPointi);

            const point& otherP = otherPatch.points()[otherPointj];
            otherPoints[pointi].append(rMergePointMap[otherPointj]);
            weights[pointi].append(1/(mag(p - otherP) + rootVSmall));
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::patchToPatchMappings::inverseDistance::generateWeights
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    generateFaceWeights
    (
        srcPatch,
        tgtPatch,
        reverse_,
        localTgtToSrc_,
        srcWeights_
    );
    generateFaceWeights
    (
        tgtPatch,
        srcPatch,
        reverse_,
        localSrcToTgt_,
        tgtWeights_
    );
}


Foam::labelList Foam::patchToPatchMappings::inverseDistance::finaliseLocal
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    generateWeights(srcPatch, tgtPatch);

    const labelList newToOld
    (
        nearby::finaliseLocal
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0
        )
    );
    tgtWeights_ = List<DynamicList<scalar>>(tgtWeights_, newToOld);

    return newToOld;
}


void Foam::patchToPatchMappings::inverseDistance::rDistributeTgt
(
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0
)
{
    // Let the base class reverse distribute the addressing
    nearby::rDistributeTgt(tgtPatch, tgtPts0);

    // Reverse distribute the face weights
    rDistributeListList
    (
        tgtPatch.size(),
        tgtMapPtr_(),
        tgtWeights_
    );
}


Foam::label Foam::patchToPatchMappings::inverseDistance::finalise
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
    const label nCouples =
        nearby::finalise
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    if (isSingleProcess())
    {
        generateWeights(srcPatch, tgtPatch);
    }

    // Normalize weights
    forAll(srcWeights_, srcFacei)
    {
        List<scalar>& ws = srcWeights_[srcFacei];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }
    forAll(tgtWeights_, tgtFacei)
    {
        List<scalar>& ws = tgtWeights_[tgtFacei];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }

    return nCouples;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::inverseDistance::inverseDistance
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool reverse
)
:
    nearby(srcPatch, tgtPatch, dict, reverse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
