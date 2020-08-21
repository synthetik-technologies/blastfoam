/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "mappedMovingPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedMovingPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedMovingPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, mappedMovingPolyPatch, dictionary);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    mappedMovingPatchBase(static_cast<const polyPatch&>(*this))
{
    //  mappedMoving is not constraint type so add mappedMoving group explicitly
    if (findIndex(inGroups(), typeName) == -1)
    {
        inGroups().append(typeName);
    }
}


Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedMovingPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm, typeName),
    mappedMovingPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedMovingPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm, typeName),
    mappedMovingPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    mappedMovingPatchBase(*this, dict)
{
    //  mappedMoving is not constraint type so add mappedMoving group explicitly
    if (findIndex(inGroups(), typeName) == -1)
    {
        inGroups().append(typeName);
    }
}


Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const mappedMovingPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    mappedMovingPatchBase(*this, pp)
{}


Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const mappedMovingPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    mappedMovingPatchBase(*this, pp)
{}


Foam::mappedMovingPolyPatch::mappedMovingPolyPatch
(
    const mappedMovingPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    mappedMovingPatchBase(*this, pp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingPolyPatch::~mappedMovingPolyPatch()
{
    mappedMovingPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedMovingPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::mappedMovingPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    polyPatch::calcGeometry(pBufs);
    mappedMovingPatchBase::clearOut();
}


void Foam::mappedMovingPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::mappedMovingPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    mappedMovingPatchBase::clearOut();
}


void Foam::mappedMovingPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::mappedMovingPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    mappedMovingPatchBase::clearOut();
}


void Foam::mappedMovingPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    mappedMovingPatchBase::write(os);
}


// ************************************************************************* //
