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
#include "triPointRef.H"
#include "patchToPatchTools.H"
#include "triangleFuncs.H"
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
    weights.setSize(otherFaces.size());
    forAll(otherFaces, facei)
    {
        weights[facei].clear();
        if (otherFaces[facei].empty()) continue;

        label otherFacei = -1;

        const labelList ovelapFaces(otherFaces[facei]);

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
                    otherPatch.localPoints()
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

        // Make sure overlapping faces are added since only one side may use
        // it for mapping

        const face& f1 = patch.localFaces()[facei];
        const vector& fc1 = patch.faceCentres()[facei];
        vector hp1, hp2;
        forAll(ovelapFaces, i)
        {
            const label otherFacej = ovelapFaces[i];
            if (findIndex(otherFaces[facei], otherFacej) >= 0)
            {
                continue;
            }

            bool hit = false;
            const face& f2 = otherPatch.localFaces()[otherFacej];
            const vector& fc2 = otherPatch.faceCentres()[otherFacej];
            forAll(f1, pi)
            {
                const point p11 = patch.localPoints()[f1[pi]];
                const point p12 = patch.localPoints()[f1[f1.fcIndex(pi)]];
                triPointRef tri1(p11, p12, fc1);

                forAll(f2, pj)
                {
                    const barycentric2D b21 = tri1.pointToBarycentric
                    (
                        otherPatch.localPoints()[f2[pj]]
                    );
                    if
                    (
                        b21.a() >= 0 && b21.a() <= 1
                     && b21.b() >= 0 && b21.b() <= 1
                     && b21.c() >= 0 && b21.c() <= 1
                    )
                    {
                        hit = true;
                        break;
                    }
                    const point p21 = tri1.barycentricToPoint(b21);
                    const point p22 = tri1.barycentricToPoint
                    (
                        tri1.pointToBarycentric
                        (
                            otherPatch.localPoints()[f2[f2.fcIndex(pj)]]
                        )
                    );

                    if
                    (
                        triangleFuncs::intersect
                        (
                            p11,
                            p12,
                            fc1,

                            p21,
                            p22,
                            fc2,

                            hp1,
                            hp2
                        )
                    )
                    {
                        hit = true;
                        break;
                    }
                }
                if (hit) break;
            }

            if (hit)
            {
                otherFaces[facei].append(otherFacej);
                weights[facei].append(0);
            }
        }
    }
}


void Foam::patchToPatchMappings::inverseDistance::generatePointWeights
(
    const primitivePatch& patch,
    const primitivePatch& otherPatch,
    const List<DynamicList<label>>& faceAddr,
    List<DynamicList<label>>& pointAddr,
    List<DynamicList<scalar>>& weights
)
{

    const pointField& points = patch.localPoints();
    const pointField& otherPoints = otherPatch.localPoints();
    const vectorField& normals = patch.faceNormals();


    labelHashSet checkedFaces;
    labelHashSet addedPoints;

    weights.setSize(points.size());
    forAll(points, pointi)
    {
        weights[pointi].clear();

        // Collect all relevant points
        const point& pt = points[pointi];
        checkedFaces.clear();
        addedPoints.clear();
        bool hit = false;

        const labelList& pointFaces = patch.pointFaces()[pointi];
        forAll(pointFaces, pfi)
        {
            const label facei = pointFaces[pfi];
            const labelList& otherFaces = faceAddr[facei];
            forAll(otherFaces, ofi)
            {
                const label otherFacei = otherFaces[ofi];
                if (checkedFaces.insert(otherFacei))
                {
                    pointHit ph = otherPatch[otherFacei].ray
                    (
                        pt,
                        normals[facei],
                        otherPoints
                    );
                    if (ph.hit())
                    {
                        hit = true;
                        const face& otherFace = otherPatch[otherFacei];
                        forAll(otherFace, pi)
                        {
                            const scalar w =
                                1.0
                               /max
                                (
                                    mag(pt - otherPoints[otherFace[pi]]),
                                    rootVSmall
                                );

                            if (!addedPoints.insert(otherFace[pi]))
                            {
                                const label curI =
                                    findIndex(pointAddr[pointi], otherFace[pi]);
                                pointAddr[pointi][curI] = otherFace[pi];
                                weights[pointi][curI] = w;
                            }
                            else
                            {
                                pointAddr[pointi].append(otherFace[pi]);
                                weights[pointi].append(w);
                            }

                        }
                    }
                    else
                    {
                        const face& otherFace = otherPatch[otherFacei];

                        forAll(otherFace, pi)
                        {
                            if (addedPoints.insert(otherFace[pi]))
                            {
                                pointAddr[pointi].append(otherFace[pi]);
                                weights[pointi].append(0.0);
                            }
                        }
                    }
                }
            }
            if (hit)
            {
                break;
            }
        }
        if (!hit)
        {
            scalar nearDistSqr = vGreat;
            label nearOtherFace = -1;
            forAll(pointAddr[pointi], opi)
            {
                const scalar distSqr =
                    magSqr(pt - otherPoints[pointAddr[pointi][opi]]);
                if (distSqr < nearDistSqr)
                {
                    nearDistSqr = distSqr;
                    nearOtherFace = opi;
                }
            }
            if (nearOtherFace >= 0)
            {
                weights[pointi][nearOtherFace] =
                    1.0/max(sqrt(nearDistSqr), small);
            }
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::patchToPatchMappings::inverseDistance::generatePointWeights
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    generatePointWeights
    (
        srcPatch,
        tgtPatch,
        localTgtFacesToSrc_,
        localTgtPointsToSrc_,
        srcPointWeights_
    );
    generatePointWeights
    (
        tgtPatch,
        srcPatch,
        localSrcFacesToTgt_,
        localSrcPointsToTgt_,
        tgtPointWeights_
    );
}


void Foam::patchToPatchMappings::inverseDistance::generateFaceWeights
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
        localTgtFacesToSrc_,
        srcFaceWeights_
    );
    generateFaceWeights
    (
        tgtPatch,
        srcPatch,
        reverse_,
        localSrcFacesToTgt_,
        tgtFaceWeights_
    );
}


Foam::labelList Foam::patchToPatchMappings::inverseDistance::finaliseLocalPoints
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    generatePointWeights(srcPatch, tgtPatch);

    const labelList newToOld
    (
        nearby::finaliseLocalPoints
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0
        )
    );
    tgtPointWeights_ = List<DynamicList<scalar>>(tgtPointWeights_, newToOld);

    return newToOld;
}



Foam::labelList Foam::patchToPatchMappings::inverseDistance::finaliseLocalFaces
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    generateFaceWeights(srcPatch, tgtPatch);

    const labelList newToOld
    (
        nearby::finaliseLocalFaces
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0
        )
    );
    tgtFaceWeights_ = List<DynamicList<scalar>>(tgtFaceWeights_, newToOld);
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
    patchToPatchTools::rDistributeListList
    (
        tgtPatch.size(),
        tgtFacesMapPtr_(),
        tgtFaceWeights_
    );

    if (needPoints_)
    {
        // Reverse distribute the point weights
        patchToPatchTools::rDistributeListList
        (
            tgtPatch.nPoints(),
            tgtPointsMapPtr_(),
            tgtPointWeights_
        );
    }
}


Foam::label Foam::patchToPatchMappings::inverseDistance::finalisePoints
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
    if (isSingleProcess())
    {
        generatePointWeights(srcPatch, tgtPatch);
    }

    const label nCouples =
        nearby::finalisePoints
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    // Remove zero weights
    DynamicList<label> newIs;
    DynamicList<scalar> newWs;
    forAll(srcPointWeights_, srcPointi)
    {
        newIs.clear();
        newWs.clear();
        labelList& Is = localTgtPointsToSrc_[srcPointi];
        scalarList& ws = srcPointWeights_[srcPointi];
        forAll(ws, i)
        {
            if (ws[i] > vSmall)
            {
                newIs.append(Is[i]);
                newWs.append(ws[i]);
            }
        }
        Is = newIs;
        ws = newWs;
    }
    forAll(tgtPointWeights_, tgtPointi)
    {
        newIs.clear();
        newWs.clear();
        labelList& Is = localSrcPointsToTgt_[tgtPointi];
        scalarList& ws = tgtPointWeights_[tgtPointi];
        forAll(ws, i)
        {
            if (ws[i] > vSmall)
            {
                newIs.append(Is[i]);
                newWs.append(ws[i]);
            }
        }
        Is = newIs;
        ws = newWs;
    }

    // Normalize weights
    forAll(srcPointWeights_, srcPointi)
    {
        scalarList& ws = srcPointWeights_[srcPointi];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }
    forAll(tgtPointWeights_, tgtPointi)
    {
        scalarList& ws = tgtPointWeights_[tgtPointi];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }

    return nCouples;
}


Foam::label Foam::patchToPatchMappings::inverseDistance::finaliseFaces
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
        nearby::finaliseFaces
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
        generateFaceWeights(srcPatch, tgtPatch);
    }

    // Normalize weights
    forAll(srcFaceWeights_, srcFacei)
    {
        List<scalar>& ws = srcFaceWeights_[srcFacei];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }
    forAll(tgtFaceWeights_, tgtFacei)
    {
        List<scalar>& ws = tgtFaceWeights_[tgtFacei];
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
    const bool needPoints,
    const bool reverse
)
:
    nearby(srcPatch, tgtPatch, dict, needPoints, reverse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
