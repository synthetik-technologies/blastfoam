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

#include "nearestPatchToPatchMapping.H"
#include "patchToPatchTools.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{
    defineTypeNameAndDebug(nearest, 0);
    addToRunTimeSelectionTable
    (
        patchToPatchMapping,
        nearest,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToPatchMappings::nearest::findNearestPoints
(
    const primitivePatch& patch,
    const primitivePatch& otherPatch,
    const List<DynamicList<label>>& faceAddr,
    List<DynamicList<label>>& pointAddr,
    List<scalar>& distances
)
{
    const pointField& points = patch.localPoints();
    const pointField& otherPoints = otherPatch.localPoints();
    const labelListList& pointFaces = patch.pointFaces();

    forAll(points, pointi)
    {
        // Collect all relevant points
        const point& pt = points[pointi];
        labelHashSet checkedPoints;

        const labelList& pFaces = pointFaces[pointi];
        forAll(pFaces, pfi)
        {
            const labelList& otherFaces = faceAddr[pFaces[pfi]];
            forAll(otherFaces, ofi)
            {
                const labelList& of = otherPatch[otherFaces[ofi]];
                forAll(of, pj)
                {
                    const label otherPointi = of[pj];
                    if (checkedPoints.insert(otherPointi))
                    {
                        pointAddr[pointi].append(otherPointi);

                        scalar distSqr =
                            magSqr(otherPoints[otherPointi] - pt);
                        if (distSqr < distances[pointi])
                        {
                            Swap
                            (
                                pointAddr[pointi].first(),
                                pointAddr[pointi].last()
                            );
                            distances[pointi] = distSqr;
                        }
                    }
                }
            }
        }
    }
}


bool Foam::patchToPatchMappings::nearest::intersectFaces
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
    if
    (
        nearby::intersectFaces
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            srcFacei,
            tgtFacei
        )
    )
    {
        const scalar dSqr =
            magSqr
            (
                srcPatch.faceCentres()[srcFacei]
              - tgtPatch.faceCentres()[tgtFacei]
            );
        if (dSqr < srcFaceDistances_[srcFacei])
        {
            srcFaceDistances_[srcFacei] = dSqr;
            Swap
            (
                localTgtFacesToSrc_[srcFacei].first(),
                localTgtFacesToSrc_[srcFacei].last()
            );
        }
        if (dSqr < tgtFaceDistances_[tgtFacei])
        {
            tgtFaceDistances_[tgtFacei] = dSqr;
            Swap
            (
                localSrcFacesToTgt_[tgtFacei].first(),
                localSrcFacesToTgt_[tgtFacei].last()
            );
        }
        return true;
    }
    return false;
}


void Foam::patchToPatchMappings::nearest::initialise
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    nearby::initialise
    (
        srcPatch,
        srcPts0,
        tgtPatch,
        tgtPts0,
        pointNormals,
        pointNormals0
    );

    srcFaceDistances_.setSize(srcPatch.size());
    srcFaceDistances_ = vGreat;

    tgtFaceDistances_.setSize(tgtPatch.size());
    tgtFaceDistances_ = vGreat;

    if (needPoints_)
    {
        srcPointDistances_.setSize(srcPatch.nPoints());
        srcPointDistances_ = vGreat;

        tgtPointDistances_.setSize(tgtPatch.nPoints());
        tgtPointDistances_ = vGreat;
    }
}


Foam::labelList Foam::patchToPatchMappings::nearest::finaliseLocalPoints
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    findNearestPoints
    (
        srcPatch,
        tgtPatch,
        localTgtFacesToSrc_,
        localTgtPointsToSrc_,
        srcPointDistances_
    );

    findNearestPoints
    (
        tgtPatch,
        srcPatch,
        localSrcFacesToTgt_,
        localSrcPointsToTgt_,
        tgtPointDistances_
    );

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
    tgtPointDistances_ = List<scalar>(tgtPointDistances_, newToOld);


    return newToOld;
}


Foam::labelList Foam::patchToPatchMappings::nearest::finaliseLocalFaces
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
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
    tgtFaceDistances_ = List<scalar>(tgtFaceDistances_, newToOld);

    return newToOld;
}


void Foam::patchToPatchMappings::nearest::rDistributeTgt
(
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0
)
{
    // Keep only the closest opposing face
    forAll(localTgtFacesToSrc_, srcFacei)
    {
        localTgtFacesToSrc_[srcFacei].resize
        (
            min(localTgtFacesToSrc_[srcFacei].size(), 1)
        );
    }
    forAll(localSrcFacesToTgt_, tgtFacei)
    {
        localSrcFacesToTgt_[tgtFacei].resize
        (
            min(localSrcFacesToTgt_[tgtFacei].size(), 1)
        );
    }

    // Create a list-list of distances to match the addressing
    List<List<scalar>> tgtFaceDistances(localSrcFacesToTgt_.size());
    forAll(localSrcFacesToTgt_, tgtFacei)
    {
        if (!localSrcFacesToTgt_[tgtFacei].empty())
        {
            tgtFaceDistances[tgtFacei].resize
            (
                1,
                tgtFaceDistances_[tgtFacei]
            );
        }
    }

    List<List<scalar>> tgtPointDistances;
    if (needPoints_)
    {
        // Keep only the closest opposing face
        forAll(localTgtPointsToSrc_, srcPointi)
        {
            localTgtPointsToSrc_[srcPointi].resize
            (
                min(localTgtPointsToSrc_[srcPointi].size(), 1)
            );
        }

        // Create a list-list of distances to match the addressing
        tgtPointDistances.setSize(localSrcPointsToTgt_.size());
        forAll(localSrcPointsToTgt_, tgtPointi)
        {
            if (!localSrcPointsToTgt_[tgtPointi].empty())
            {
                localSrcPointsToTgt_[tgtPointi].resize(1);
                tgtPointDistances[tgtPointi].resize
                (
                    1,
                    tgtPointDistances_[tgtPointi]
                );
            }
        }
    }

    // Let the base class reverse distribute the addressing
    nearby::rDistributeTgt(tgtPatch, tgtPts0);


    // Reverse distribute the face distances
    patchToPatchTools::rDistributeListList
    (
        tgtPatch.size(),
        tgtFacesMapPtr_(),
        tgtFaceDistances
    );

    // If there is more than one address, remove all but the closest
    tgtFaceDistances_.resize(localSrcFacesToTgt_.size());
    forAll(localSrcFacesToTgt_, tgtFacei)
    {
        if (localSrcFacesToTgt_[tgtFacei].size() > 1)
        {
            const label neari = findMin(tgtFaceDistances[tgtFacei]);
            const label srcFacei = localSrcFacesToTgt_[tgtFacei][neari];

            localSrcFacesToTgt_[tgtFacei].setSize(1);
            localSrcFacesToTgt_[tgtFacei][0] = srcFacei;
            tgtFaceDistances_[tgtFacei] = tgtFaceDistances[tgtFacei][neari];
        }
    }

    if (needPoints_)
    {
        // Reverse distribute the face distances
        patchToPatchTools::rDistributeListList
        (
            tgtPatch.nPoints(),
            tgtPointsMapPtr_(),
            tgtPointDistances
        );

        // If there is more than one address, remove all but the closest
        tgtPointDistances_.resize(localSrcPointsToTgt_.size());
        forAll(localSrcPointsToTgt_, tgtPointi)
        {
            if (localSrcPointsToTgt_[tgtPointi].size() > 1)
            {
                const label neari = findMin(tgtPointDistances[tgtPointi]);
                const label srcPointi = localSrcPointsToTgt_[tgtPointi][neari];

                localSrcPointsToTgt_[tgtPointi].resize(1);
                localSrcPointsToTgt_[tgtPointi][0] = srcPointi;
                tgtPointDistances_[tgtPointi] = tgtPointDistances[tgtPointi][neari];
            }
        }
    }
}


Foam::label Foam::patchToPatchMappings::nearest::finalisePoints
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
        findNearestPoints
        (
            srcPatch,
            tgtPatch,
            localTgtFacesToSrc_,
            localTgtPointsToSrc_,
            srcPointDistances_
        );
        findNearestPoints
        (
            tgtPatch,
            srcPatch,
            localSrcFacesToTgt_,
            localSrcPointsToTgt_,
            tgtPointDistances_
        );
    }

    // Keep only the closest opposing point
    srcPointWeights_.setSize(localTgtPointsToSrc_.size());
    forAll(localTgtPointsToSrc_, srcPointi)
    {
        localTgtPointsToSrc_[srcPointi].resize
        (
            min(localTgtPointsToSrc_[srcPointi].size(), 1)
        );
        srcPointWeights_[srcPointi].resize
        (
            localTgtPointsToSrc_[srcPointi].size(),
            1.0
        );
    }
    tgtPointWeights_.setSize(localSrcPointsToTgt_.size());
    forAll(localSrcPointsToTgt_, tgtPointi)
    {
        localSrcPointsToTgt_[tgtPointi].resize
        (
            min(localSrcPointsToTgt_[tgtPointi].size(), 1)
        );
        tgtPointWeights_[tgtPointi].resize
        (
            localSrcPointsToTgt_[tgtPointi].size(),
            1.0
        );
    }

    return
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
}


Foam::label Foam::patchToPatchMappings::nearest::finaliseFaces
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
    // Keep only the closest opposing face
    srcFaceWeights_.setSize(localTgtFacesToSrc_.size());
    forAll(localTgtFacesToSrc_, srcFacei)
    {
        localTgtFacesToSrc_[srcFacei].resize
        (
            min(localTgtFacesToSrc_[srcFacei].size(), 1)
        );
        srcFaceWeights_[srcFacei].resize
        (
            localTgtFacesToSrc_[srcFacei].size(),
            1.0
        );
    }

    tgtFaceWeights_.setSize(localSrcFacesToTgt_.size());
    forAll(localSrcFacesToTgt_, tgtFacei)
    {
        localSrcFacesToTgt_[tgtFacei].resize
        (
            min(localSrcFacesToTgt_[tgtFacei].size(), 1)
        );
        tgtFaceWeights_[tgtFacei].resize
        (
            localSrcFacesToTgt_[tgtFacei].size(),
            1.0
        );
    }

    return
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::nearest::nearest
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


Foam::patchToPatchMappings::nearest::nearest
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool needPoints,
    const bool reverse,
    const scalar minBbDim
)
:
    nearby(srcPatch, tgtPatch, dict, needPoints, reverse, minBbDim)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
