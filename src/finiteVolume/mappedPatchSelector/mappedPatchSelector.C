/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License

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

#include "mappedPatchBase.H"

#include "mappedMovingPatchBase.H"
#include "mappedMovingWallPolyPatch.H"
#include "mappedMovingVariableThicknessWallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchSelector, 0);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchSelector::mappedPatchSelector
(
    const polyPatch& pp
)
:
    patch_(pp),
    mappedPatchPtr_(nullptr),
    mappedMovingPatchPtr_(nullptr)
{
    if (isA<mappedPatchBase>(pp))
    {
        mappedPatchPtr_ = &refCast<const mappedPatchBase>(pp);
    }
    else if (isA<mappedMovingPatchBase>(pp))
    {
        mappedMovingPatchPtr_ =
            &refCast<const mappedMovingPatchBase>(pp);
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


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::mappedPatchSelector::clearOut()
{
    if (mappedPatchPtr_)
    {
        //- Clear this patch
        mappedPatchBase& mpb =
            const_cast<mappedPatchBase&>(*mappedPatchPtr_);
        mpb.clearOut();

        // Clear the sample patch
        polyPatch& spp(const_cast<polyPatch&>(samplePolyPatch()));
        dynamicCast<mappedPatchBase>(spp).clearOut();
    }

    if (mappedMovingPatchPtr_)
    {
        //- Clear this patch
        mappedMovingPatchBase& mmpb =
            const_cast<mappedMovingPatchBase&>(*mappedMovingPatchPtr_);
        mmpb.clearOut();

        // Clear the sample patch
        polyPatch& spp(const_cast<polyPatch&>(samplePolyPatch()));
        dynamicCast<mappedMovingPatchBase>(spp).clearOut();
    }

    if (pointInterpolatorPtr_.valid())
    {
        pointInterpolatorPtr_->movePoints();
    }
}

// ************************************************************************* //
