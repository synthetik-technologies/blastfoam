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

#include "burstProcessorCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstProcessorCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, burstProcessorCyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstProcessorCyclicPolyPatch::burstProcessorCyclicPolyPatch
(
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const word& referPatchName,
    const word& patchType
)
:
    processorCyclicPolyPatch
    (
        size,
        start,
        index,
        bm,
        myProcNo,
        neighbProcNo,
        referPatchName,
        patchType
    )
{}


Foam::burstProcessorCyclicPolyPatch::burstProcessorCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    processorCyclicPolyPatch(name, dict, index, bm, patchType)
{}


Foam::burstProcessorCyclicPolyPatch::burstProcessorCyclicPolyPatch
(
    const burstProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    processorCyclicPolyPatch(pp, bm)
{}


Foam::burstProcessorCyclicPolyPatch::burstProcessorCyclicPolyPatch
(
    const burstProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    processorCyclicPolyPatch(pp, bm, index, newSize, newStart)
{}


Foam::burstProcessorCyclicPolyPatch::burstProcessorCyclicPolyPatch
(
    const burstProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& referPatchName
)
:
    processorCyclicPolyPatch(pp, bm, index, newSize, newStart, referPatchName)
{}


Foam::burstProcessorCyclicPolyPatch::burstProcessorCyclicPolyPatch
(
    const burstProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    processorCyclicPolyPatch(pp, bm, index, mapAddressing, newStart)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstProcessorCyclicPolyPatch::~burstProcessorCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
