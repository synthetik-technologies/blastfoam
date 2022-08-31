/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "oversetPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, oversetPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, oversetPolyPatch, dictionary);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::oversetPolyPatch::oversetPolyPatch
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
    masterPatchID_(-1)
{
    //  'overset' is not constraint type so add to group explicitly
    if (findIndex(inGroups(), typeName) < 0)
    {
        inGroups().append(typeName);
    }
}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    masterPatchID_(-1)
{
    //  'overset' is not constraint type so add to group explicitly
    if (findIndex(inGroups(), typeName) < 0)
    {
        inGroups().append(typeName);
    }
}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const oversetPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    masterPatchID_(-1)
{}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const oversetPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    masterPatchID_(-1)
{}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const oversetPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    masterPatchID_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetPolyPatch::~oversetPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::oversetPolyPatch::master() const
{
    if (masterPatchID_ == -1)
    {
        // Find lowest numbered overset patch
        const polyBoundaryMesh& pbm = boundaryMesh();

        forAll(pbm, patchi)
        {
            if (isA<oversetPolyPatch>(pbm[patchi]))
            {
                masterPatchID_ = patchi;
                break;
            }
        }

        if (masterPatchID_ != -1 && masterPatchID_ > 0)
        {
            WarningInFunction<< "The master overset patch is not the"
                << " first patch. Generally the first patch should be an"
                << " overset patch to guarantee consistent operation." << endl;
        }
    }

    return index() == masterPatchID_;
}


// ************************************************************************* //
