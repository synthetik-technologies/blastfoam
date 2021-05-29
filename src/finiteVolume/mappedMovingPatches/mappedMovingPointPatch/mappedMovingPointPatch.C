/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "mappedMovingPointPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedMovingPointPatch, 0);

    addToRunTimeSelectionTable(facePointPatch, mappedMovingPointPatch, polyPatch);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingPointPatch::mappedMovingPointPatch
(
    const polyPatch& pp,
    const pointBoundaryMesh& bm
)
:
    facePointPatch(pp, bm),
    mappedMovingPointPatchBase(pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingPointPatch::~mappedMovingPointPatch()
{
    mappedMovingPointPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedMovingPointPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    facePointPatch::initCalcGeometry(pBufs);
}


void Foam::mappedMovingPointPatch::calcGeometry(PstreamBuffers& pBufs)
{
    facePointPatch::calcGeometry(pBufs);
    mappedMovingPointPatchBase::clearOut();
}


void Foam::mappedMovingPointPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    facePointPatch::initMovePoints(pBufs, p);
}


void Foam::mappedMovingPointPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    facePointPatch::movePoints(pBufs, p);
    mappedMovingPointPatchBase::clearOut();
}


void Foam::mappedMovingPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
}


void Foam::mappedMovingPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    mappedMovingPointPatchBase::clearOut();
}


void Foam::mappedMovingPointPatch::write(Ostream& os) const
{}


// ************************************************************************* //
