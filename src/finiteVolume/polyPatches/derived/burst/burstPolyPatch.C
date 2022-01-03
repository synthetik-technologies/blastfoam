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
    polyPatch::initCalcGeometry(pBufs);
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
    polyPatch::calcGeometry(pBufs);
}


void Foam::burstPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::burstPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
}


void Foam::burstPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    needMap_ = true;
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::burstPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    needMap_ = true;
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
    polyPatch(name, size, start, index, bm, patchType),
    intact_(size, 1.0),
    usePressure_(true),
    pBurst_(great),
    pRef_(size, 0.0),
    useImpulse_(false),
    impulseBurst_(great),
    partialBurst_(false),
    needMap_(false)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::burstPolyPatch::burstPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    intact_(this->size(), 1.0),
    usePressure_(dict.lookupOrDefault("usePressure", true)),
    pBurst_(great),
    pRef_(this->size(), 0.0),
    useImpulse_(dict.lookupOrDefault("useImpulse", false)),
    impulseBurst_(great),
    partialBurst_(dict.lookupOrDefault("partialBurst", false)),
    needMap_(false)
{
    if (usePressure_)
    {
        pBurst_ = dict.lookupOrDefault<scalar>("pBurst", 0.0);
    }
    if (useImpulse_)
    {
        impulseBurst_ = dict.lookup<scalar>("impulseBurst");
    }
}


Foam::burstPolyPatch::burstPolyPatch
(
    const burstPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    intact_(pp.intact_),
    usePressure_(pp.usePressure_),
    pBurst_(pp.pBurst_),
    pRef_(pp.pRef_),
    useImpulse_(pp.useImpulse_),
    impulseBurst_(pp.impulseBurst_),
    partialBurst_(pp.partialBurst_),
    needMap_(false)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::burstPolyPatch::burstPolyPatch
(
    const burstPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    intact_(pp.intact_),
    usePressure_(pp.usePressure_),
    pBurst_(pp.pBurst_),
    pRef_(pp.pRef_),
    useImpulse_(pp.useImpulse_),
    impulseBurst_(pp.impulseBurst_),
    partialBurst_(pp.partialBurst_),
    needMap_(false)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::burstPolyPatch::burstPolyPatch
(
    const burstPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    intact_(pp.intact_),
    usePressure_(pp.usePressure_),
    pBurst_(pp.pBurst_),
    pRef_(pp.pRef_),
    useImpulse_(pp.useImpulse_),
    impulseBurst_(pp.impulseBurst_),
    partialBurst_(pp.partialBurst_),
    needMap_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPolyPatch::~burstPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::burstPolyPatch::needMap() const
{
    return needMap_;
}


void Foam::burstPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
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
