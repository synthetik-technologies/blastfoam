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

#include "burstCyclicAMIPolyPatch.H"
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
    defineTypeNameAndDebug(burstCyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, burstCyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, burstCyclicAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::burstCyclicAMIPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::initCalcGeometry(pBufs);
}


void Foam::burstCyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::calcGeometry(pBufs);
}


void Foam::burstCyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);
}


void Foam::burstCyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicAMIPolyPatch::movePoints(pBufs, p);
}


void Foam::burstCyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    needMap_ = true;
    cyclicAMIPolyPatch::initUpdateMesh(pBufs);
}


void Foam::burstCyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::updateMesh(pBufs);
    needMap_ = true;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstCyclicAMIPolyPatch::burstCyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const bool AMIRequireMatch,
    const AMIInterpolation::interpolationMethod AMIMethod
)
:
    cyclicAMIPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        AMIRequireMatch,
        AMIMethod
    ),
    intact_(size, 1.0),
    usePressure_(true),
    pBurst_(great),
    useImpulse_(false),
    impulseBurst_(great),
    partialBurst_(false),
    needMap_(false)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::burstCyclicAMIPolyPatch::burstCyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const bool AMIRequireMatch,
    const AMIInterpolation::interpolationMethod AMIMethod
)
:
    cyclicAMIPolyPatch
    (
        name,
        dict,
        index,
        bm,
        patchType,
        AMIRequireMatch,
        AMIMethod
    ),
    intact_(this->size(), 1.0),
    usePressure_(dict.lookupOrDefault("usePressure", true)),
    pBurst_(great),
    useImpulse_(dict.lookupOrDefault("useImpulse", false)),
    impulseBurst_(great),
    partialBurst_(dict.lookupOrDefault("partialBurst", false)),
    needMap_(false)
{
    if (usePressure_)
    {
        pBurst_ = dict.lookup<scalar>("pBurst");
    }
    if (useImpulse_)
    {
        impulseBurst_ = dict.lookup<scalar>("impulseBurst");
    }
}


Foam::burstCyclicAMIPolyPatch::burstCyclicAMIPolyPatch
(
    const burstCyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    intact_(pp.intact_),
    usePressure_(pp.usePressure_),
    pBurst_(pp.pBurst_),
    useImpulse_(pp.useImpulse_),
    impulseBurst_(pp.impulseBurst_),
    partialBurst_(pp.partialBurst_),
    needMap_(false)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::burstCyclicAMIPolyPatch::burstCyclicAMIPolyPatch
(
    const burstCyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    intact_(pp.intact_),
    usePressure_(pp.usePressure_),
    pBurst_(pp.pBurst_),
    useImpulse_(pp.useImpulse_),
    impulseBurst_(pp.impulseBurst_),
    partialBurst_(pp.partialBurst_),
    needMap_(false)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::burstCyclicAMIPolyPatch::burstCyclicAMIPolyPatch
(
    const burstCyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    intact_(pp.intact_),
    usePressure_(pp.usePressure_),
    pBurst_(pp.pBurst_),
    useImpulse_(pp.useImpulse_),
    impulseBurst_(pp.impulseBurst_),
    partialBurst_(pp.partialBurst_),
    needMap_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstCyclicAMIPolyPatch::~burstCyclicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::burstCyclicAMIPolyPatch::needMap() const
{
    return needMap_;
}


void Foam::burstCyclicAMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);
    writeEntry(os, "usePressure", usePressure_);
    if (usePressure_)
    {
        writeEntry(os, "pBurst", pBurst_);
    }
    writeEntry(os, "useImpulse", useImpulse_);
    if (useImpulse_)
    {
        writeEntry(os, "impulseBurst", impulseBurst_);
    }
    writeEntry(os, "partialBurst", partialBurst_);
}


// ************************************************************************* //
