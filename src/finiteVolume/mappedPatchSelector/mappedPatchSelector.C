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

#include "mappedPatchSelector.H"
#include "mappedWallFvPatch.H"
#include "mappedMovingWallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchSelector, 0);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchSelector::mappedPatchSelector
(
    const fvPatch& pp
)
:
    mappedPatchPtr_(nullptr),
    mappedMovingPatchPtr_(nullptr)
{
    if (isA<mappedWallFvPatch>(pp))
    {
        mappedPatchPtr_ = &refCast<const mappedPatchBase>(pp.patch());
    }
    else
    {
        mappedMovingPatchPtr_ = &refCast<const mappedMovingPatchBase>(pp.patch());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchSelector::~mappedPatchSelector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedPatchSelector::sampleMesh() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->sampleMesh();
    }

    return mappedMovingPatchPtr_->sampleMesh();
}


const Foam::polyPatch& Foam::mappedPatchSelector::samplePolyPatch() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->samplePolyPatch();
    }

    return mappedMovingPatchPtr_->samplePolyPatch();
}


Foam::tmp<Foam::pointField> Foam::mappedPatchSelector::samplePoints() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->samplePoints();
    }

    return mappedMovingPatchPtr_->samplePoints();
}


Foam::pointIndexHit Foam::mappedPatchSelector::facePoint
(
    const polyMesh& mesh,
    const label facei,
    const polyMesh::cellDecomposition decompMode
)
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->facePoint(mesh, facei, decompMode);
    }

    return mappedMovingPatchPtr_->facePoint(mesh, facei, decompMode);
}


// ************************************************************************* //
