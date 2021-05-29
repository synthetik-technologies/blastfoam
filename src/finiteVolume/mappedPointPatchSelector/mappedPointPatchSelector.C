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

#include "mappedPointPatchSelector.H"

#include "mappedPointPatch.H"

#include "mappedMovingPointPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPointPatchSelector, 0);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPointPatchSelector::mappedPointPatchSelector
(
    const pointPatch& pp
)
:
    mappedMovingPatchPtr_(nullptr)
{
    if (isA<mappedMovingPointPatch>(pp))
    {
        mappedMovingPatchPtr_ =
            &refCast<const mappedMovingPointPatchBase>(pp);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPointPatchSelector::~mappedPointPatchSelector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedPointPatchSelector::sampleMesh() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->sampleMesh();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->sampleMesh();
}


const Foam::polyPatch& Foam::mappedPointPatchSelector::samplePolyPatch() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->samplePolyPatch();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->samplePolyPatch();
}


Foam::tmp<Foam::pointField> Foam::mappedPointPatchSelector::samplePoints() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->samplePoints();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->samplePoints();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //
// Function to clear patches if an update has occurred
void Foam::mappedPointPatchSelector::clearMappedPatches(fvMesh& mesh)
{
    pointMesh& pMesh = mesh.lookupObjectRef<pointMesh>("pointMesh");
    forAll(pMesh.boundary(), patchi)
    {
        if (isA<mappedMovingPointPatch>(pMesh.boundary()[patchi]))
        {
            pointBoundaryMesh& pbMesh =
                const_cast<pointBoundaryMesh&>(pMesh.boundary());
            refCast<mappedMovingPointPatchBase>
            (
                pbMesh[patchi]
            ).clearOut();
        }
    }
}


bool Foam::mappedPointPatchSelector::setMappedPatchDisplacement
(
    fvMesh& ownMesh,
    const fvMesh& neiMesh
)
{
    pointMesh& pMesh = ownMesh.lookupObjectRef<pointMesh>("pointMesh");
    if (neiMesh.foundObject<pointVectorField>("pointD"))
    {
        const pointVectorField& pointD =
            neiMesh.lookupObject<pointVectorField>("pointD");
        forAll(pMesh.boundary(), patchi)
        {
            if (isA<mappedMovingPointPatch>(pMesh.boundary()[patchi]))
            {
                pointBoundaryMesh& pbMesh =
                    const_cast<pointBoundaryMesh&>(pMesh.boundary());
                mappedMovingPointPatchBase& mmppb =
                    refCast<mappedMovingPointPatchBase>
                    (
                        pbMesh[patchi]
                    );
                if (mmppb.sampleRegion() == neiMesh.name())
                {
                    mmppb.setOffsets(pointD);
                }
            }
        }
        return true;
    }
    return false;
}


bool Foam::mappedPointPatchSelector::isAMappedType
(
    const pointPatch& patch
)
{
    if
    (
        isA<mappedMovingPointPatch>(patch)
     || isA<mappedPointPatch>(patch)
    )
    {
        return true;
    }
    return false;
}


void Foam::mappedPointPatchSelector::replaceMappedWithMappedMoving
(
    fvMesh& mesh
)
{
    pointMesh& pMesh = mesh.lookupObjectRef<pointMesh>("pointMesh");
    pointBoundaryMesh& bm =
        const_cast<pointBoundaryMesh&>(pMesh.boundary());

    forAll(bm, patchi)
    {
        if (isA<mappedPointPatch>(bm[patchi]))
        {
            bm.set
            (
                patchi,
                new mappedMovingPointPatch
                (
                    mesh.boundaryMesh()[patchi],
                    bm
                )
            );
        }
    }

    mesh.clearOut();
}


// ************************************************************************* //
