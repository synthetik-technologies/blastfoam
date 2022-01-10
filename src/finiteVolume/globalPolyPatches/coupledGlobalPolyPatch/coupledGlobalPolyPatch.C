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
    if (patchToPatchInterpPtr_)
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    const pointField& pts = this->globalPatch().points();
    const pointField& samplePts = samplePatch().globalPatch().points();
    boundBox bb(pts, false);
    boundBox sampleBb(samplePts, false);

    boundBox inflatedBb(bb);
    boundBox inflatedSampleBb(sampleBb);
    inflatedBb.inflate(1e-4);
    inflatedSampleBb.inflate(1e-4);

    if (!inflatedBb.contains(sampleBb))
    {
        labelHashSet unmapped;
        forAll(samplePts, pti)
        {
            if (!inflatedBb.contains(samplePts[pti]))
            {
                unmapped.insert(pti);
            }
        }
        samplePatch().unmappedPoints_ = unmapped.toc();

        unmapped.clear();
        const pointField& fc = samplePatch().globalPatch().faceCentres();
        forAll(fc, fi)
        {
            if (!inflatedBb.contains(fc[fi]))
            {
                unmapped.insert(fi);
            }
        }
        samplePatch().unmappedFaces_ = unmapped.toc();
    }
    if (!inflatedSampleBb.contains(bb))
    {
        labelHashSet unmapped;
        forAll(pts, pti)
        {
            if (!inflatedSampleBb.contains(pts[pti]))
            {
                unmapped.insert(pti);
            }
        }
        unmappedPoints_ = unmapped.toc();

        unmapped.clear();
        const pointField& fc = this->globalPatch().faceCentres();
        forAll(fc, fi)
        {
            if (!inflatedSampleBb.contains(fc[fi]))
            {
                unmapped.insert(fi);
            }
        }
        unmappedFaces_ = unmapped.toc();
    }

    patchToPatchInterpPtr_ =
        patchToPatchMapping::New
        (
            dict_.found("mappingType")
          ? dict_
          : samplePatch().dict_,
            patch(),
            samplePatch().patch(),
            *this,
            samplePatch()
        ).ptr();
    samplePatch().setPatchToPatchInterp(patchToPatchInterpPtr_);
}


void Foam::coupledGlobalPolyPatch::setPatchToPatchInterp
(
    patchToPatchMapping* interpPtr
) const
{
    patchToPatchInterpPtr_ = interpPtr;
}


void Foam::coupledGlobalPolyPatch::clearOut() const
{
    clearInterp();
    globalPolyPatch::clearOut();
}


void Foam::coupledGlobalPolyPatch::clearInterp(const bool top) const
{
    if (patchToPatchInterpPtr_)
    {
        deleteDemandDrivenData
        (
            patchToPatchInterpPtr_
        );
        unmappedFaces_.clear();
        unmappedPoints_.clear();

        samplePatch().setPatchToPatchInterp(nullptr);
        samplePatch().unmappedFaces_.clear();
        samplePatch().unmappedPoints_.clear();
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
    dict_(dict),
    sampleRegion_(dict.lookup("sampleRegion")),
    samplePatch_(dict.lookup("samplePatch")),
    patchToPatchInterpPtr_(nullptr)
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

const Foam::patchToPatchMapping&
Foam::coupledGlobalPolyPatch::patchToPatchInterpolator() const
{
    if (!patchToPatchInterpPtr_)
    {
        calcPatchToPatchInterp();
    }
    return *patchToPatchInterpPtr_;
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
