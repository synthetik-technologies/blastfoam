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

#include "burstProcessorCyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "edgeList.H"
#include "globalPolyBoundaryMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstProcessorCyclicPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        burstProcessorCyclicPointPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::burstProcessorCyclicPointPatch::initCalcGeometry(PstreamBuffers&)
{}


void Foam::burstProcessorCyclicPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::burstProcessorCyclicPointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void Foam::burstProcessorCyclicPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::burstProcessorCyclicPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    burstProcessorCyclicPointPatch::initCalcGeometry(pBufs);
    curTimeIndex_ = -1;
}


void Foam::burstProcessorCyclicPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    burstProcessorCyclicPointPatch::calcGeometry(pBufs);
    curTimeIndex_ = -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstProcessorCyclicPointPatch::burstProcessorCyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    processorCyclicPointPatch(patch, bm),
    burstProcessorCyclicPolyPatch_(refCast<const burstProcessorCyclicPolyPatch>(patch)),
    intact_(patch.nPoints(), 1.0),
    curTimeIndex_(-1)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstProcessorCyclicPointPatch::~burstProcessorCyclicPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstProcessorCyclicPointPatch::update(const label curTimeIndex)
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
                burstProcessorCyclicPolyPatch_.boundaryMesh().mesh().thisDb()
            )
        )[burstProcessorCyclicPolyPatch_].faceToPoint(burstProcessorCyclicPolyPatch_.intact());
}


// ************************************************************************* //
