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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledGlobalPolyPatch, 0);
    defineRunTimeSelectionTable(coupledGlobalPolyPatch, patch);
    addToRunTimeSelectionTable
    (
        coupledGlobalPolyPatch,
        coupledGlobalPolyPatch,
        patch
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::patchToPatchMapping&
Foam::coupledGlobalPolyPatch::patchToPatchInterpolator
(
    const bool needPoints
) const
{
    if (needPoints && !needPoints_)
    {
        needPoints_ = true;
        clearInterp();
    }

    if (!patchToPatchInterpPtr_)
    {
        calcPatchToPatchInterp();
    }
    return *patchToPatchInterpPtr_;
}

void Foam::coupledGlobalPolyPatch::calcPatchToPatchInterp() const
{
    if (patchToPatchInterpPtr_)
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    isSrc_ = dict_.found("mappingType");
    samplePatch().isSrc_ = !isSrc_;

    const coupledGlobalPolyPatch& masterPatch =
        isSrc_
      ? *this
      : samplePatch();
    const coupledGlobalPolyPatch& slavePatch =
        isSrc_
      ? samplePatch()
      : *this;

    patchToPatchInterpPtr_ =
        patchToPatchMapping::New
        (
            masterPatch.physicalPatch(),
            slavePatch.physicalPatch(),
            masterPatch.dict_,
            needPoints_,
            false
        ).ptr();

    bool masterMoving =
        &masterPatch.physicalPatch() != &masterPatch.physicalPatch0();
    bool slaveMoving =
        &slavePatch.physicalPatch() != &slavePatch.physicalPatch0();

    patchToPatchInterpPtr_->update
    (
        masterMoving
      ? masterPatch.physicalPatch0().points()
      : NullObjectRef<pointField>(),
        slaveMoving
      ? slavePatch.physicalPatch0().points()
      : NullObjectRef<pointField>(),
        masterPatch.physicalPatch().pointNormals(),
        masterMoving
      ? masterPatch.physicalPatch0().pointNormals()
      : NullObjectRef<vectorField>()
    );

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
        deleteDemandDrivenData(patchToPatchInterpPtr_);
        deleteDemandDrivenData(unmappedFacesPtr_);
        deleteDemandDrivenData(unmappedPointsPtr_);

        samplePatch().setPatchToPatchInterp(nullptr);
        deleteDemandDrivenData(samplePatch().unmappedFacesPtr_);
        deleteDemandDrivenData(samplePatch().unmappedPointsPtr_);
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
    needPoints_(false),
    sampleRegion_(dict.lookup("sampleRegion")),
    samplePatch_(dict.lookup("samplePatch")),
    patchToPatchInterpPtr_(nullptr),
    unmappedFacesPtr_(nullptr),
    unmappedPointsPtr_(nullptr)
{}


Foam::autoPtr<Foam::coupledGlobalPolyPatch> Foam::coupledGlobalPolyPatch::New
(
    const dictionary& dict,
    const polyPatch& patch
)
{
    patchConstructorTable::iterator cstrIter =
        patchConstructorTablePtr_->find(patch.type());

    if (cstrIter != patchConstructorTablePtr_->end())
    {
        return cstrIter()(dict, patch);
    }
    return autoPtr<coupledGlobalPolyPatch>
    (
        new coupledGlobalPolyPatch(dict, patch)
    );
}

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


const Foam::labelList&
Foam::coupledGlobalPolyPatch::unmappedFaces() const
{
    if (!unmappedFacesPtr_)
    {
        unmappedFacesPtr_ = new labelList
        (
            patchToPatchInterpolator().unmappedFaces(physicalPatch())
        );
    }
    return *unmappedFacesPtr_;
}


const Foam::labelList&
Foam::coupledGlobalPolyPatch::unmappedPoints() const
{
    if (!unmappedPointsPtr_)
    {
        unmappedPointsPtr_ = new labelList
        (
            patchToPatchInterpolator().unmappedPoints(physicalPatch())
        );
    }
    return *unmappedPointsPtr_;
}


void Foam::coupledGlobalPolyPatch::update()
{
    globalPolyPatch::update();
}


void Foam::coupledGlobalPolyPatch::movePoints(const bool clear)
{
    globalPolyPatch::movePoints(clear);
}


void Foam::coupledGlobalPolyPatch::updateMesh()
{
    clearOut();
}


// ************************************************************************* //
