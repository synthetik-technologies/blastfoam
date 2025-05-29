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

#include "nearbyPatchToPatchMapping.H"
#include "boundSphere.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchMappings
{
    defineTypeNameAndDebug(nearby, 0);
//     addToRunTimeSelectionTable(patchToPatchMapping, dictionary, nearby);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::nearby::nearby
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool needPoints,
    const bool reverse
)
:
    nearby
    (
        srcPatch,
        tgtPatch,
        dict,
        needPoints,
        reverse,
        dict.lookupOrDefault("minBbDim", -1.0)
    )
{}


Foam::patchToPatchMappings::nearby::nearby
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool needPoints,
    const bool reverse,
    const scalar minBbDim
)
:
    patchToPatchMapping(srcPatch, tgtPatch, dict, needPoints, reverse),
    minBbDim_(minBbDim),
    srcSpheres_(0),
    tgtSpheres_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::nearby::~nearby()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatchMappings::nearby::makeBb
(
    const face& f,
    const pointField& pts,
    const vectorField& pns
) const
{
    const treeBoundBox bb(pts, f);

    const point c = bb.midpoint();
    const scalar l = minBbDim_ < 0 ? bb.maxDim() : max(minBbDim_, bb.maxDim());

    return treeBoundBox(c - l*vector::one, c + l*vector::one);
}


bool Foam::patchToPatchMappings::nearby::intersectFaces
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
    const point& centreA = srcSpheres_[srcFacei].first();
    const scalar& radiusA = srcSpheres_[srcFacei].second();
    const point& centreB = tgtSpheres_[tgtFacei].first();
    const scalar& radiusB = tgtSpheres_[tgtFacei].second();

    if (magSqr(centreA - centreB) < sqr(radiusA + radiusB))
    {
        localTgtFacesToSrc_[srcFacei].append(tgtFacei);
        localSrcFacesToTgt_[tgtFacei].append(srcFacei);
        return true;
    }
    return false;
}


void Foam::patchToPatchMappings::nearby::initialise
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    patchToPatchMapping::initialise
    (
        srcPatch,
        srcPts0,
        tgtPatch,
        tgtPts0,
        pointNormals,
        pointNormals0
    );

    srcSpheres_.setSize(srcPatch.size());
    forAll(srcPatch, srcFacei)
    {
        srcSpheres_[srcFacei] =
            boundSphere
            (
                UIndirectList<point>(srcPatch.localPoints(), srcPatch[srcFacei])
            );
        srcSpheres_[srcFacei].second() = max(srcSpheres_[srcFacei].second(), minBbDim_);
    }

    tgtSpheres_.setSize(tgtPatch.size());
    forAll(tgtPatch, tgtFacei)
    {
        tgtSpheres_[tgtFacei] =
            boundSphere
            (
                UIndirectList<point>(tgtPatch.localPoints(), tgtPatch[tgtFacei])
            );
        tgtSpheres_[tgtFacei].second() = max(tgtSpheres_[tgtFacei].second(), minBbDim_);
    }
}

// ************************************************************************* //
