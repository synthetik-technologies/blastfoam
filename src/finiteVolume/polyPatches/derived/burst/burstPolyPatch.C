/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "burstPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "matchPoints.H"
#include "EdgeMap.H"
#include "Time.H"
#include "transformField.H"
#include "SubField.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, burstPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, burstPolyPatch, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::burstPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::initCalcGeometry(pBufs);
}


void Foam::burstPolyPatch::initCalcGeometry
(
    const primitivePatch& referPatch,
    pointField& nbrCtrs,
    vectorField& nbrAreas,
    pointField& nbrCc
)
{}


void Foam::burstPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::calcGeometry(pBufs);
}


void Foam::burstPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::initMovePoints(pBufs, p);
}


void Foam::burstPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::movePoints(pBufs, p);
}


void Foam::burstPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::initUpdateMesh(pBufs);
}


void Foam::burstPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::updateMesh(pBufs);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstPolyPatch::burstPolyPatch
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
    burstPolyPatchBase(dynamicCast<const polyPatch>(*this))
{}


Foam::burstPolyPatch::burstPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    burstPolyPatchBase(*this, dict)
{}


Foam::burstPolyPatch::burstPolyPatch
(
    const burstPolyPatch& bpp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(bpp, bm),
    burstPolyPatchBase(*this, bpp)
{}


Foam::burstPolyPatch::burstPolyPatch
(
    const burstPolyPatch& bpp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(bpp, bm, index, newSize, newStart),
    burstPolyPatchBase(*this, bpp, newSize, newStart)
{}


Foam::burstPolyPatch::burstPolyPatch
(
    const burstPolyPatch& bpp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    wallPolyPatch(bpp, bm, index, mapAddressing, newStart),
    burstPolyPatchBase(*this, bpp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPolyPatch::~burstPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    burstPolyPatchBase::write(os);
}


// ************************************************************************* //
