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
        if (dSqr < srcDistances_[srcFacei])
        {
            srcDistances_[srcFacei] = dSqr;
            Swap
            (
                localTgtToSrc_[srcFacei].first(),
                localTgtToSrc_[srcFacei].last()
            );
        }
        if (dSqr < tgtDistances_[tgtFacei])
        {
            tgtDistances_[tgtFacei] = dSqr;
            Swap
            (
                localSrcToTgt_[tgtFacei].first(),
                localSrcToTgt_[tgtFacei].last()
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

    srcDistances_.setSize(srcPatch.size());
    srcDistances_ = vGreat;

    tgtDistances_.setSize(tgtPatch.size());
    tgtDistances_ = vGreat;

    forAll(srcWeights_, i)
    {
        srcWeights_[i].setSize(1);
        srcWeights_[i][0] = 1.0;
    }

    forAll(tgtWeights_, i)
    {
        tgtWeights_[i].setSize(1);
        tgtWeights_[i][0] = 1.0;
    }
}


Foam::labelList Foam::patchToPatchMappings::nearest::finaliseLocal
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
    tgtDistances_ = List<scalar>(tgtDistances_, newToOld);

    return newToOld;
}


void Foam::patchToPatchMappings::nearest::rDistributeTgt
(
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0
)
{
    // Keep only the closest opposing face
    forAll(localTgtToSrc_, srcFacei)
    {
        localTgtToSrc_[srcFacei].resize
        (
            min(localTgtToSrc_[srcFacei].size(), 1)
        );
    }
    forAll(localSrcToTgt_, tgtFacei)
    {
        localSrcToTgt_[tgtFacei].resize
        (
            min(localSrcToTgt_[tgtFacei].size(), 1)
        );
    }

    // Create a list-list of distances to match the addressing
    List<List<scalar>> tgtDistances(localSrcToTgt_.size());
    forAll(localSrcToTgt_, tgtFacei)
    {
        if (!localSrcToTgt_[tgtFacei].empty())
        {
            tgtDistances[tgtFacei].resize(1, tgtDistances_[tgtFacei]);
        }
    }


    // Let the base class reverse distribute the addressing
    nearby::rDistributeTgt(tgtPatch, tgtPts0);

    // Reverse distribute the face distances
    rDistributeListList
    (
        tgtPatch.size(),
        tgtMapPtr_(),
        tgtDistances
    );

    // If there is more than one address, remove all but the closest
    tgtDistances_.resize(localSrcToTgt_.size());
    forAll(localSrcToTgt_, tgtFacei)
    {
        if (localSrcToTgt_[tgtFacei].size() > 1)
        {
            const label neari = findMin(tgtDistances[tgtFacei]);

            localSrcToTgt_[tgtFacei].resize(1);
            localSrcToTgt_[tgtFacei][0] = localSrcToTgt_[tgtFacei][neari];
            tgtDistances_[tgtFacei] = tgtDistances[tgtFacei][neari];
        }
    }
}


Foam::label Foam::patchToPatchMappings::nearest::finalise
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
    forAll(localTgtToSrc_, srcFacei)
    {
        localTgtToSrc_[srcFacei].resize
        (
            min(localTgtToSrc_[srcFacei].size(), 1)
        );
    }
    forAll(localSrcToTgt_, tgtFacei)
    {
        localSrcToTgt_[tgtFacei].resize
        (
            min(localSrcToTgt_[tgtFacei].size(), 1)
        );
    }

    return
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::nearest::nearest
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
