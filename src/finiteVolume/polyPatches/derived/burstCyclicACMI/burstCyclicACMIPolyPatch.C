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

#include "burstCyclicACMIPolyPatch.H"
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
    defineTypeNameAndDebug(burstCyclicACMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, burstCyclicACMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, burstCyclicACMIPolyPatch, dictionary);
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstCyclicACMIPolyPatch::burstCyclicACMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicACMIPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType
    ),
    burstPolyPatchBase(dynamicCast<const polyPatch>(*this))
{}


Foam::burstCyclicACMIPolyPatch::burstCyclicACMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicACMIPolyPatch
    (
        name,
        dict,
        index,
        bm,
        patchType
    ),
    burstPolyPatchBase(*this, dict)
{}


Foam::burstCyclicACMIPolyPatch::burstCyclicACMIPolyPatch
(
    const burstCyclicACMIPolyPatch& bcpp,
    const polyBoundaryMesh& bm
)
:
    cyclicACMIPolyPatch(bcpp, bm),
    burstPolyPatchBase(*this, bcpp)
{}


Foam::burstCyclicACMIPolyPatch::burstCyclicACMIPolyPatch
(
    const burstCyclicACMIPolyPatch& bcpp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& nonOverlapPatchName
)
:
    cyclicACMIPolyPatch
    (
        bcpp,
        bm,
        index,
        newSize,
        newStart,
        nbrPatchName,
        nonOverlapPatchName
    ),
    burstPolyPatchBase(*this, bcpp, newSize, newStart)
{}


Foam::burstCyclicACMIPolyPatch::burstCyclicACMIPolyPatch
(
    const burstCyclicACMIPolyPatch& bcpp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicACMIPolyPatch(bcpp, bm, index, mapAddressing, newStart),
    burstPolyPatchBase(*this, bcpp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstCyclicACMIPolyPatch::~burstCyclicACMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstCyclicACMIPolyPatch::write(Ostream& os) const
{
    cyclicACMIPolyPatch::write(os);
    burstPolyPatchBase::write(os);
}


// ************************************************************************* //
