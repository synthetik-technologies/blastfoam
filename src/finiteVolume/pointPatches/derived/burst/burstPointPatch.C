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

#include "burstPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "edgeList.H"
#include "globalPolyBoundaryMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        burstPointPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::burstPointPatch::initCalcGeometry(PstreamBuffers&)
{}


void Foam::burstPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::burstPointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void Foam::burstPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::burstPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    burstPointPatch::initCalcGeometry(pBufs);
    curTimeIndex_ = -1;
}


void Foam::burstPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    burstPointPatch::calcGeometry(pBufs);
    curTimeIndex_ = -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstPointPatch::burstPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    facePointPatch(patch, bm),
    burstPolyPatch_(refCast<const burstPolyPatch>(patch)),
    intact_(patch.nPoints(), 1.0),
    curTimeIndex_(-1)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPointPatch::~burstPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstPointPatch::update(const label curTimeIndex)
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
                burstPolyPatch_.boundaryMesh().mesh().thisDb()
            )
        )[burstPolyPatch_].faceToPoint(burstPolyPatch_.intact());
}


// ************************************************************************* //
