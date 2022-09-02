/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "cellBox.H"
#include "mapDistribute.H"

namespace Foam
{
namespace processorLODs
{
    defineTypeNameAndDebug(cellBox, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::boundBox Foam::processorLODs::cellBox::calcSrcBox
(
    const label srcObji
) const
{
    const UList<label>& cellFaces = srcCells_[srcObji];

    boundBox bb(srcPoints_, srcFaces_[cellFaces[0]], false);
    for (label i = 1; i < cellFaces.size(); ++i)
    {
        boundBox tmpBb(srcPoints_, srcFaces_[cellFaces[i]]);
        bb.min() = min(bb.min(), tmpBb.min());
        bb.max() = min(bb.max(), tmpBb.max());
    }

    return bb;
}


Foam::boundBox Foam::processorLODs::cellBox::calcTgtBox
(
    const label tgtObji
) const
{
    const UList<label>& cellFaces = tgtCells_[tgtObji];

    boundBox bb(tgtPoints_, tgtFaces_[cellFaces[0]], false);
    for (label i = 1; i < cellFaces.size(); ++i)
    {
        boundBox tmpBb(tgtPoints_, tgtFaces_[cellFaces[i]]);
        bb.min() = min(bb.min(), tmpBb.min());
        bb.max() = min(bb.max(), tmpBb.max());
    }

    return bb;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::processorLODs::cellBox::cellBox
(
    const cellList& srcCells,
    const faceList& srcFaces,
    const UList<point>& srcPoints,
    const cellList& tgtCells,
    const faceList& tgtFaces,
    const UList<point>& tgtPoints,
    const label maxObjectsPerLeaf,
    const label nObjectsOfType,
    const label nRefineIterMax
)
:
    faceBox
    (
        srcFaces,
        srcPoints,
        tgtFaces,
        tgtPoints,
        maxObjectsPerLeaf,
        nObjectsOfType,
        nRefineIterMax
    ),
    srcCells_(srcCells),
    tgtCells_(tgtCells)
{}


Foam::autoPtr<Foam::mapDistribute> Foam::processorLODs::cellBox::map()
{
    return createMap(srcCells_.size(), tgtCells_.size());
}


// ************************************************************************* //
