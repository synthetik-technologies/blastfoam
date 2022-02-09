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

#include "burstCyclicPolyPatch.H"
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
    defineTypeNameAndDebug(burstCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, burstCyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, burstCyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

// void Foam::burstCyclicPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
// {
//     polyPatch::initCalcGeometry(pBufs);
// }
//
//
// void Foam::burstCyclicPolyPatch::initCalcGeometry
// (
//     const primitivePatch& referPatch,
//     pointField& nbrCtrs,
//     vectorField& nbrAreas,
//     pointField& nbrCc
// )
// {
//     cyclicPolyPatch::initCalcGeometry(referPatch, nbrCtrs, nbrAreas, nbrCc);
// }
//
//
// void Foam::burstCyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
// {
//     cyclicPolyPatch::calcGeometry(pBufs);
// //     static_cast<cyclicTransform&>(*this) =
// //         cyclicTransform(true);
// }
//
//
// void Foam::burstCyclicPolyPatch::initMovePoints
// (
//     PstreamBuffers& pBufs,
//     const pointField& p
// )
// {
//     polyPatch::initMovePoints(pBufs, p);
// }
//
//
// void Foam::burstCyclicPolyPatch::movePoints
// (
//     PstreamBuffers& pBufs,
//     const pointField& p
// )
// {
//     polyPatch::movePoints(pBufs, p);
// }
//
//
// void Foam::burstCyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
// {
//     cyclicPolyPatch::initUpdateMesh(pBufs);
// }
//
//
// void Foam::burstCyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
// {
//     cyclicPolyPatch::updateMesh(pBufs);
//     burstPolyPatchBase::updateMesh();
// }


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstCyclicPolyPatch::burstCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, size, start, index, bm, patchType),
    burstPolyPatchBase(static_cast<const polyPatch&>(*this))
{}


Foam::burstCyclicPolyPatch::burstCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& nbrPatchName
)
:
    cyclicPolyPatch(name, size, start, index, bm, patchType, nbrPatchName),
    burstPolyPatchBase(static_cast<const polyPatch&>(*this))
{}


Foam::burstCyclicPolyPatch::burstCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, dict, index, bm, patchType),
    burstPolyPatchBase(*this, dict)
{
    static_cast<cyclicTransform&>(*this) = cyclicTransform(dict, true);
}


Foam::burstCyclicPolyPatch::burstCyclicPolyPatch
(
    const burstCyclicPolyPatch& bcpp,
    const polyBoundaryMesh& bm
)
:
    cyclicPolyPatch(bcpp, bm),
    burstPolyPatchBase(*this, bcpp)
{}


Foam::burstCyclicPolyPatch::burstCyclicPolyPatch
(
    const burstCyclicPolyPatch& bcpp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& neiName
)
:
    cyclicPolyPatch(bcpp, bm, index, newSize, newStart, neiName),
    burstPolyPatchBase(*this, bcpp, newSize, newStart)
{}


Foam::burstCyclicPolyPatch::burstCyclicPolyPatch
(
    const burstCyclicPolyPatch& bcpp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicPolyPatch(bcpp, bm, index, mapAddressing, newStart),
    burstPolyPatchBase(*this, bcpp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstCyclicPolyPatch::~burstCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstCyclicPolyPatch::write(Ostream& os) const
{
    cyclicPolyPatch::write(os);
    burstPolyPatchBase::write(os);
}


// ************************************************************************* //
