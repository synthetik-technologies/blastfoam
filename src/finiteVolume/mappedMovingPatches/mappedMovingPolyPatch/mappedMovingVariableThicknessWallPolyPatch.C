/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "mappedMovingVariableThicknessWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedMovingVariableThicknessWallPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedMovingVariableThicknessWallPolyPatch,
        word
    );

    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedMovingVariableThicknessWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingVariableThicknessWallPolyPatch::mappedMovingVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    mappedMovingWallPolyPatch(name, size, start, index, bm, patchType),
    thickness_(size)
{}


Foam::mappedMovingVariableThicknessWallPolyPatch::mappedMovingVariableThicknessWallPolyPatch
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
    mappedMovingWallPolyPatch
    (
        name,
        size,
        start,
        index,
        sampleRegion,
        samplePatch,
        bm
    ),
    thickness_(size)
{}


Foam::mappedMovingVariableThicknessWallPolyPatch::mappedMovingVariableThicknessWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    mappedMovingWallPolyPatch(name, dict, index, bm, patchType),
    thickness_(scalarField("thickness", dict, this->size()))
{}


Foam::mappedMovingVariableThicknessWallPolyPatch::
mappedMovingVariableThicknessWallPolyPatch
(
    const mappedMovingVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    mappedMovingWallPolyPatch(pp, bm),
    thickness_(pp.thickness_)
{}


Foam::mappedMovingVariableThicknessWallPolyPatch::mappedMovingVariableThicknessWallPolyPatch
(
    const mappedMovingVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    mappedMovingWallPolyPatch(pp, bm, index, newSize, newStart),
    thickness_(newSize)
{}


Foam::mappedMovingVariableThicknessWallPolyPatch::mappedMovingVariableThicknessWallPolyPatch
(
    const mappedMovingVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    mappedMovingWallPolyPatch(pp, bm, index, mapAddressing, newStart),
    thickness_(pp.size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingVariableThicknessWallPolyPatch::
~mappedMovingVariableThicknessWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedMovingVariableThicknessWallPolyPatch::
write(Foam::Ostream& os) const
{
    writeEntry(os, "thickness", thickness_);
}


// ************************************************************************* //
