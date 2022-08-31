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

#include "faceBox.H"
#include "mapDistribute.H"

namespace Foam
{
namespace processorLODs
{
    defineTypeNameAndDebug(faceBox, 0);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::boundBox Foam::processorLODs::faceBox::calcSrcBox
(
    const label srcObji
) const
{
    return boundBox(srcPoints_, srcFaces_[srcObji], false);
}


Foam::boundBox Foam::processorLODs::faceBox::calcTgtBox
(
    const label tgtObji
) const
{
    return boundBox(tgtPoints_, tgtFaces_[tgtObji], false);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::processorLODs::faceBox::faceBox
(
    const faceList& srcFaces,
    const UList<point>& srcPoints,
    const faceList& tgtFaces,
    const UList<point>& tgtPoints,
    const label maxObjectsPerLeaf,
    const label nObjectsOfType,
    const label nRefineIterMax
)
:
    box(srcPoints, tgtPoints, maxObjectsPerLeaf, nObjectsOfType),
    srcFaces_(srcFaces),
    tgtFaces_(tgtFaces)
{}


Foam::autoPtr<Foam::mapDistribute> Foam::processorLODs::faceBox::map()
{
    return createMap(srcFaces_.size(), tgtFaces_.size());
}


// ************************************************************************* //
