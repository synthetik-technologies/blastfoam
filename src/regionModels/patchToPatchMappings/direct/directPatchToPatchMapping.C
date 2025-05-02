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

#include "directPatchToPatchMapping.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchMappings
{

    defineTypeNameAndDebug(direct, 0);
    addToRunTimeSelectionTable(patchToPatchMapping, direct, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::patchToPatchMappings::direct::forwardCheck
(
    const word& elemType,
    const vectorField& pts,
    const List<DynamicList<label>>& localOther,
    const bool isSrc
)
{
    forAll(localOther, i)
    {
        if (localOther[i].size() != 1)
        {
            FatalErrorInFunction
                << (isSrc ? "Source " : "Target ")
                << elemType << " #" << i << " at "
                << pts[i]
                << " did not match a " << elemType << " on the "
                << (isSrc ? "target" : "source")
                << " side" << exit(FatalError);
            return false;
        }
    }
    return true;
};

// Make sure every face is referenced by exactly one face
bool Foam::patchToPatchMappings::direct::reverseCheck
(
    const word& elemType,
    const vectorField& pts,
    const List<DynamicList<label>>& otherLocal,
    const autoPtr<distributionMap>& mapPtr,
    const bool isSrc
)
{
    labelList count
    (
        mapPtr.valid() ? mapPtr->constructSize() : pts.size(),
        0
    );

    forAll(otherLocal, i)
    {
        forAll(otherLocal[i], j)
        {
            count[otherLocal[i][j]]++;
        }
    }

    if (mapPtr.valid())
    {
        distributionMapBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            pts.size(),
            mapPtr->constructMap(),
            false,
            mapPtr->subMap(),
            false,
            count,
            plusEqOp<label>(),
            flipOp(),
            label(0)
        );
    }

    forAll(count, i)
    {
        if (count[i] != 1)
        {
            FatalErrorInFunction
                << (isSrc ? "Source " : "Target ")
                << elemType << " #" << i << " at "
                << pts[i]
                << " did not match a " << elemType << " on the "
                << (isSrc ? "target" : "source")
                << " side" << exit(FatalError);

            return false;
        }
    }
    return true;
};

Foam::label Foam::patchToPatchMappings::direct::finalisePoints
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
        nearest::finalisePoints
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    forwardCheck("point", srcPatch.localPoints(), localTgtPointsToSrc_, true);
    forwardCheck("point", tgtPatch.localPoints(), localSrcPointsToTgt_, false);

    reverseCheck("point", srcPatch.localPoints(), localSrcPointsToTgt_, srcPointsMapPtr_, true);
    reverseCheck("point", tgtPatch.localPoints(), localTgtPointsToSrc_, tgtPointsMapPtr_, false);

    return nCouples;
}

Foam::label Foam::patchToPatchMappings::direct::finaliseFaces
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
        nearest::finaliseFaces
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    forwardCheck("face", srcPatch.faceCentres(), localTgtFacesToSrc_, true);
    forwardCheck("face", tgtPatch.faceCentres(), localSrcFacesToTgt_, false);

    reverseCheck("face", srcPatch.faceCentres(), localSrcFacesToTgt_, srcFacesMapPtr_, true);
    reverseCheck("face", tgtPatch.faceCentres(), localTgtFacesToSrc_, tgtFacesMapPtr_, false);

    return nCouples;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::direct::direct
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool needPoints,
    const bool reverse
)
:
    nearest(srcPatch, tgtPatch, dict, needPoints, reverse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
