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

#include "mappedMovingWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedMovingPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedMovingWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedMovingWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedMovingWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingWallPolyPatch::mappedMovingWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, size, start, index, bm, patchType),
    mappedMovingPatchBase(static_cast<const polyPatch&>(*this))
{
    //  mappedMoving is not constraint type so add mappedMoving group explicitly
    if (findIndex(inGroups(), mappedMovingPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedMovingPolyPatch::typeName);
    }
}


Foam::mappedMovingWallPolyPatch::mappedMovingWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const word& samplePatch,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm, typeName),
    mappedMovingPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        samplePatch
    )
{}


Foam::mappedMovingWallPolyPatch::mappedMovingWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    mappedMovingPatchBase(*this, dict)
{
    //  mappedMoving is not constraint type so add mappedMoving group explicitly
    if (findIndex(inGroups(), mappedMovingPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedMovingPolyPatch::typeName);
    }
}


Foam::mappedMovingWallPolyPatch::mappedMovingWallPolyPatch
(
    const mappedMovingWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    mappedMovingPatchBase(*this, pp)
{}


Foam::mappedMovingWallPolyPatch::mappedMovingWallPolyPatch
(
    const mappedMovingWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    mappedMovingPatchBase(*this, pp)
{}


Foam::mappedMovingWallPolyPatch::mappedMovingWallPolyPatch
(
    const mappedMovingWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, mapAddressing, newStart),
    mappedMovingPatchBase(*this, pp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingWallPolyPatch::~mappedMovingWallPolyPatch()
{
    mappedMovingPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedMovingWallPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::initCalcGeometry(pBufs);
}


void Foam::mappedMovingWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::calcGeometry(pBufs);
    mappedMovingPatchBase::clearOut();
}


void Foam::mappedMovingWallPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::initMovePoints(pBufs, p);
}


void Foam::mappedMovingWallPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::movePoints(pBufs, p);
    mappedMovingPatchBase::clearOut();
}


void Foam::mappedMovingWallPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::initUpdateMesh(pBufs);
}


void Foam::mappedMovingWallPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::updateMesh(pBufs);
    mappedMovingPatchBase::clearOut();
}


void Foam::mappedMovingWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    mappedMovingPatchBase::write(os);
}


// ************************************************************************* //
