/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cohesivePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "patchZones.H"
#include "matchPoints.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cohesivePolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, cohesivePolyPatch, word);
addToRunTimeSelectionTable(polyPatch, cohesivePolyPatch, dictionary);


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

cohesivePolyPatch::cohesivePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType)
{}


cohesivePolyPatch::cohesivePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType)
{}


cohesivePolyPatch::cohesivePolyPatch
(
    const cohesivePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm)
{}


cohesivePolyPatch::cohesivePolyPatch
(
    const cohesivePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart)
{}


Foam::cohesivePolyPatch::cohesivePolyPatch
(
    const cohesivePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cohesivePolyPatch::~cohesivePolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cohesivePolyPatch::initOrder
(
    PstreamBuffers&,
    const primitivePatch& pp
) const
{}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool cohesivePolyPatch::order
(
    PstreamBuffers&,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return false;


    // Grab faceMap
    SortableList<label> topoFaceMap(faceMap);

    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (pp.size() == 0)
    {
        // No faces, nothing to change.
        return false;
    }

    label oldFacesStartIndex = -1;
    for (label i = 0; i < topoFaceMap.size(); i += 2)
    {
        if (topoFaceMap[i + 1] - topoFaceMap[i] == 1)
        {
            oldFacesStartIndex = i;
            break;
        }
    }

    label oldFacesSize = 0;
    if (oldFacesStartIndex != -1)
    {
        oldFacesSize = oldFacesSize + 2;
        for (label i = oldFacesStartIndex + 2; i < topoFaceMap.size(); i += 2)
        {
            oldFacesSize += 2*(topoFaceMap[i + 1]-topoFaceMap[i]);
        }
    }

    // Indices of faces on half0
    labelList half0ToPatch(pp.size());
    // Indices of faces on half1
    labelList half1ToPatch(pp.size());

    {
        // Calculate normals
        vectorField normals(pp.size());

        forAll (pp, faceI)
        {
            normals[faceI] = pp[faceI].normal(pp.points());
        }

        normals /= mag(normals) + VSMALL;

        label n0Faces = 0;
        label n1Faces = 0;

        label sizeByTwo = oldFacesSize/2;

        if (oldFacesStartIndex != -1)
        {
            for
            (
                label i = oldFacesStartIndex;
                i < oldFacesStartIndex+sizeByTwo;
                i++
            )
            {
                half0ToPatch[n0Faces++] = i;
                half1ToPatch[n1Faces++] = i + sizeByTwo;
            }

            for (label i = 0; i < oldFacesStartIndex; i += 2)
            {
                half0ToPatch[n0Faces++] = i;
                half1ToPatch[n1Faces++] = i + 1;
            }

            for
            (
                label i=oldFacesStartIndex+oldFacesSize;
                i<topoFaceMap.size();
                i = i + 2
            )
            {
                half0ToPatch[n0Faces++] = i;
                half1ToPatch[n1Faces++] = i + 1;
            }
        }
        else
        {
            for
            (
                label i = 0;
                i < topoFaceMap.size();
                i += 2
            )
            {
                half0ToPatch[n0Faces++] = i;
                half1ToPatch[n1Faces++] = i + 1;
            }
        }

        half0ToPatch.setSize(n0Faces);
        half1ToPatch.setSize(n1Faces);

        Pout<< "cohesivePolyPatch::order : "
            << "Number of faces per zone:("
            << n0Faces << ' ' << n1Faces << ')' << endl;
    }

    if (half0ToPatch.size() != half1ToPatch.size())
    {
        SeriousErrorIn
        (
            "cohesivePolyPatch::order"
            "(const primitivePatch&, labelList&, labelList&) const"
        )   << " patch:" << name() << " : "
            << "Patch " << name() << " gets decomposed in two zones of"
            << "inequal size: " << half0ToPatch.size()
            << " and " << half1ToPatch.size() << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }

    forAll (half0ToPatch, faceI)
    {
        // Label in original patch
        label patchFaceI = topoFaceMap.indices()[half0ToPatch[faceI]];

        faceMap[patchFaceI] = faceI;
    }

    forAll (half1ToPatch, faceI)
    {
        // Label in original patch
        label patchFaceI = topoFaceMap.indices()[half1ToPatch[faceI]];

        faceMap[patchFaceI] = half0ToPatch.size() + faceI;
    }

    forAll (faceMap, faceI)
    {
        if (faceMap[faceI] != faceI)
        {
            return true;
        }
    }

    return false;
}


void cohesivePolyPatch::syncOrder() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
