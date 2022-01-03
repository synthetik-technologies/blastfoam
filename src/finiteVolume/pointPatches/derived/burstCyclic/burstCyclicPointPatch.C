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

#include "burstCyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "edgeList.H"
#include "globalPolyBoundaryMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstCyclicPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        burstCyclicPointPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::burstCyclicPointPatch::initCalcGeometry(PstreamBuffers&)
{}


void Foam::burstCyclicPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::burstCyclicPointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void Foam::burstCyclicPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::burstCyclicPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    burstCyclicPointPatch::initCalcGeometry(pBufs);
    curTimeIndex_ = -1;
}


void Foam::burstCyclicPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    burstCyclicPointPatch::calcGeometry(pBufs);
    curTimeIndex_ = -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstCyclicPointPatch::burstCyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    cyclicPointPatch(patch, bm),
    burstCyclicPolyPatch_(refCast<const burstCyclicPolyPatch>(patch)),
    intact_(patch.nPoints(), 1.0),
    curTimeIndex_(-1)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstCyclicPointPatch::~burstCyclicPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstCyclicPointPatch::update(const label curTimeIndex)
{
    if (curTimeIndex_ == curTimeIndex)
    {
        return;
    }
    curTimeIndex_ = curTimeIndex;

    intact_ =
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                burstCyclicPolyPatch_.boundaryMesh().mesh().thisDb()
            )
        )[burstCyclicPolyPatch_].faceToPoint(burstCyclicPolyPatch_.intact());
}


// ************************************************************************* //
