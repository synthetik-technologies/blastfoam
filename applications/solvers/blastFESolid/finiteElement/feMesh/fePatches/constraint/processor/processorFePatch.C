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

#include "processorFePatch.H"
#include "feBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "feMesh.H"
#include "faceList.H"
#include "primitiveFacePatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorFePatch, 0);
    addToRunTimeSelectionTable
    (
        fePatch,
        processorFePatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::processorFePatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    if (boundaryMesh().mesh().order() < 2)
    {
        // Algorithm:
        // Depending on whether the patch is a master or a slave, get the primitive
        // patch points and filter away the points from the global patch.

        // Create the reversed patch and pick up its points
        // so that the order is correct
        const polyPatch& pp = this->patch();

        faceList masterFaces(pp.size());

        forAll(pp, facei)
        {
            masterFaces[facei] = pp[facei].reverseFace();
        }

        reverseMeshNodes_ = primitiveFacePatch
        (
            masterFaces,
            pp.points()
        ).meshPoints();
    }
}


void Foam::processorFePatch::calcGeometry(PstreamBuffers& pBufs)
{}


void Foam::processorFePatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::processorFePatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::processorFePatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    fePatch::initUpdateMesh(pBufs);
    processorFePatch::initCalcGeometry(pBufs);
}


void Foam::processorFePatch::updateMesh(PstreamBuffers& pBufs)
{
    fePatch::updateMesh(pBufs);
    processorFePatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorFePatch::processorFePatch
(
    const polyPatch& patch,
    const feBoundaryMesh& bm
)
:
    coupledFePatch(patch, bm),
    procPolyPatch_(refCast<const processorPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorFePatch::~processorFePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::processorFePatch::reverseMeshNodes() const
{
    return reverseMeshNodes_;
}


// ************************************************************************* //
