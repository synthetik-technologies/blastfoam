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
    const primitivePatch& patch,
    const List<DynamicList<label>>& localOtherFaces,
    const bool isSrc
)
{
    forAll(localOtherFaces, facei)
    {
        if (localOtherFaces[facei].size() != 1)
        {
            FatalErrorInFunction
                << (isSrc ? "Source" : "Target")
                << " face #" << facei << " at "
                << patch.faceCentres()[facei]
                << " did not match a face on the "
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
    const primitivePatch& patch,
    const List<DynamicList<label>>& otherLocalFaces,
    const autoPtr<mapDistribute>& mapPtr,
    const bool isSrc
)
{
    labelList count
    (
        mapPtr.valid() ? mapPtr->constructSize() : patch.size(),
        0
    );

    forAll(otherLocalFaces, otherFacei)
    {
        forAll(otherLocalFaces[otherFacei], i)
        {
            count[otherLocalFaces[otherFacei][i]] ++;
        }
    }

    if (mapPtr.valid())
    {
        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            patch.size(),
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

    forAll(count, facei)
    {
        if (count[facei] != 1)
        {
            FatalErrorInFunction
                << (isSrc ? "Source" : "Target")
                << " face #" << facei << " at "
                << patch.faceCentres()[facei]
                << " did not match a face on the "
                << (isSrc ? "target" : "source")
                << " side" << exit(FatalError);

            return false;
        }
    }
    return true;
};

Foam::label Foam::patchToPatchMappings::direct::finalise
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
        nearest::finalise
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    forwardCheck(srcPatch, localTgtToSrc_, true);
    forwardCheck(tgtPatch, localSrcToTgt_, false);

    reverseCheck(srcPatch, localSrcToTgt_, srcMapPtr_, true);
    reverseCheck(tgtPatch, localTgtToSrc_, tgtMapPtr_, false);

    return nCouples;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::direct::direct
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool reverse
)
:
    nearest(srcPatch, tgtPatch, dict, reverse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
