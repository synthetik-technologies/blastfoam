/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coupledGlobalPolyPatch.H"
#include "globalPolyBoundaryMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledGlobalPolyPatch, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coupledGlobalPolyPatch::calcPatchToPatchInterp() const
{
    if (patchToSampleInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    patchToSampleInterpPtr_.set
    (
        new standAlonePatchToPatchInterpolation
        (
            globalPatch(),
            samplePatch().globalPatch()/*,
            intersection::algorithm::visible,
            intersection::direction::contactSphere*/
        )
    );
}


void Foam::coupledGlobalPolyPatch::clearOut() const
{
    clearInterp();
    globalPolyPatch::clearOut();
}


void Foam::coupledGlobalPolyPatch::clearInterp(const bool top) const
{
    if (patchToSampleInterpPtr_.valid())
    {
        patchToSampleInterpPtr_.clear();
    }
    if (top)
    {
        samplePatch().clearInterp(false);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::coupledGlobalPolyPatch::coupledGlobalPolyPatch
(
    const dictionary& dict,
    const polyPatch& patch
)
:
    globalPolyPatch(dict, patch),
    sampleRegion_(dict.lookup("sampleRegion")),
    samplePatch_(dict.lookup("samplePatch")),
    patchToSampleInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledGlobalPolyPatch::~coupledGlobalPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::polyMesh& Foam::coupledGlobalPolyPatch::sampleMesh() const
{
    return mesh_.time().lookupObject<polyMesh>(sampleRegion_);
}


const Foam::coupledGlobalPolyPatch&
Foam::coupledGlobalPolyPatch::samplePatch() const
{
    return globalPolyBoundaryMesh::New(sampleMesh())(samplePatch_);
}

const Foam::coupledGlobalPolyPatch::standAlonePatchToPatchInterpolation&
Foam::coupledGlobalPolyPatch::patchToPatchInterpolator() const
{
    if (!patchToSampleInterpPtr_.valid())
    {
        calcPatchToPatchInterp();
    }
    return patchToSampleInterpPtr_();
}


void Foam::coupledGlobalPolyPatch::movePoints()
{
    clearOut();
}


void Foam::coupledGlobalPolyPatch::updateMesh()
{
    clearOut();
}


// ************************************************************************* //
