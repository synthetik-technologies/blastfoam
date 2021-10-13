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
#include "volFields.H"
#include "pointFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledGlobalPolyPatch, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::standAlonePatch& Foam::coupledGlobalPolyPatch::dispPatch() const
{
    if (dispPatchPtr_.valid())
    {
        return dispPatchPtr_();
    }

    // Create globalPatch with optional displacement
    if (this->mesh_.foundObject<volVectorField>(displacementField_))
    {
        dispPatchPtr_.set
        (
            new standAlonePatch
            (
                this->globalPatch(),
                this->globalPatch().localPoints()
              - this->interpolator().faceToPointInterpolate
                (
                    this->patchFaceToGlobal
                    (
                        this->mesh_.lookupObject<volVectorField>
                        (
                            displacementField_
                        ).boundaryField()[this->patch_.index()]
                    )
                )
            )
        );
    }
    else if (this->mesh_.foundObject<pointVectorField>(displacementField_))
    {
        dispPatchPtr_.set
        (
            new standAlonePatch
            (
                this->globalPatch(),
                this->globalPatch().localPoints()
              - this->patchFaceToGlobal
                (
                    this->mesh_.lookupObject<pointVectorField>
                    (
                        displacementField_
                    ).boundaryField()
                    [this->patch_.index()].patchInternalField()
                )
            )
        );
    }
    else
    {
        return globalPatch();
    }
    return dispPatchPtr_();
}


const Foam::standAlonePatch&
Foam::coupledGlobalPolyPatch::dispSamplePatch() const
{
    if (dispSamplePatchPtr_.valid())
    {
        return dispPatchPtr_();
    }

    //- Construct optional displaced global patch for the sample patch
    const polyMesh& nbrMesh = sampleMesh();
    const globalPolyPatch& nbrPatch =
        globalPolyBoundaryMesh::New(nbrMesh)[samplePatch_];

    if (nbrMesh.foundObject<volVectorField>(sampleDisplacementField_))
    {
        dispSamplePatchPtr_.set
        (
            new standAlonePatch
            (
                nbrPatch.globalPatch(),
                nbrPatch.interpolator().faceToPointInterpolate
                (
                    nbrPatch.patchFaceToGlobal
                    (
                        nbrMesh.lookupObject<volVectorField>
                        (
                            sampleDisplacementField_
                        ).boundaryField()[nbrPatch.patch().index()]
                    )
                ) + nbrPatch.globalPatch().localPoints()
            )
        );
    }
    else if (nbrMesh.foundObject<pointVectorField>(sampleDisplacementField_))
    {
        dispSamplePatchPtr_.set
        (
            new standAlonePatch
            (
                nbrPatch.globalPatch(),
                nbrPatch.patchPointToGlobal
                (
                    nbrMesh.lookupObject<pointVectorField>
                    (
                        sampleDisplacementField_
                    ).boundaryField()
                    [nbrPatch.patch().index()].patchInternalField()
                ) + nbrPatch.globalPatch().localPoints()
            )
        );
    }
    else
    {
        return nbrPatch.globalPatch();
    }
    return dispSamplePatchPtr_();
}


void Foam::coupledGlobalPolyPatch::calcPatchToPatchInterp() const
{
    if (patchToSampleInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    patchToSampleInterpPtr_.reset
    (
        new standAlonePatchToPatchInterpolation
        (
            dispPatch(),
            dispSamplePatch()
        )
    );
}


void Foam::coupledGlobalPolyPatch::clearOut() const
{
    globalPolyPatch::clearOut();
    clearInterp();
}


void Foam::coupledGlobalPolyPatch::clearInterp() const
{
    if
    (
        dispPatchPtr_.valid()
     || dispSamplePatchPtr_.valid()
     || patchToSampleInterpPtr_.valid()
    )
    {
        dispPatchPtr_.clear();
        dispSamplePatchPtr_.clear();
        patchToSampleInterpPtr_.clear();
        samplePatch().clearInterp();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::coupledGlobalPolyPatch::coupledGlobalPolyPatch
(
    const word& patchName,
    const polyMesh& mesh,
    const word& displacementField,
    const word& sampleDisplacementField
)
:
    globalPolyPatch(patchName, mesh),
    sampleRegion_
    (
        dynamicCast<const mappedPatchBase>(patch_).sampleRegion()
    ),
    samplePatch_
    (
        dynamicCast<const mappedPatchBase>(patch_).sampleRegion()
    ),
    displacementField_(displacementField),
    sampleDisplacementField_(sampleDisplacementField),
    dispPatchPtr_(nullptr),
    dispSamplePatchPtr_(nullptr),
    patchToSampleInterpPtr_(nullptr)
{}


// Construct from dictionary
Foam::coupledGlobalPolyPatch::coupledGlobalPolyPatch
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    globalPolyPatch(dict, mesh),
    sampleRegion_
    (
        dynamicCast<const mappedPatchBase>(patch_).sampleRegion()
    ),
    samplePatch_
    (
        dynamicCast<const mappedPatchBase>(patch_).sampleRegion()
    ),
    displacementField_
    (
        dict.lookupOrDefault("displacementField", word::null)
    ),
    sampleDisplacementField_
    (
        dict.lookupOrDefault("sampleDisplacementField", word::null)
    ),
    dispPatchPtr_(nullptr),
    dispSamplePatchPtr_(nullptr),
    patchToSampleInterpPtr_(nullptr)
{}


Foam::coupledGlobalPolyPatch::coupledGlobalPolyPatch
(
    const polyPatch& pp,
    const word& displacementField,
    const word& sampleDisplacementField
)
:
    globalPolyPatch(pp),
    sampleRegion_
    (
        dynamicCast<const mappedPatchBase>(patch_).sampleRegion()
    ),
    samplePatch_
    (
        dynamicCast<const mappedPatchBase>(patch_).samplePatch()
    ),
    displacementField_(displacementField),
    sampleDisplacementField_(sampleDisplacementField),
    dispPatchPtr_(nullptr),
    dispSamplePatchPtr_(nullptr),
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
    return dynamicCast<const coupledGlobalPolyPatch>
    (
        globalPolyBoundaryMesh::New(sampleMesh())[samplePatch_]
    );
}

const Foam::coupledGlobalPolyPatch::standAlonePatchToPatchInterpolation&
Foam::coupledGlobalPolyPatch::patchToPatchInterpolator() const
{
    patchToSampleInterpPtr_.clear();
    if (!patchToSampleInterpPtr_.valid())
    {
        calcPatchToPatchInterp();
    }
    return patchToSampleInterpPtr_();
}


void Foam::coupledGlobalPolyPatch::movePoints()
{
    globalPolyPatch::clearOut();
    clearOut();
}


void Foam::coupledGlobalPolyPatch::updateMesh()
{
    globalPolyPatch::clearOut();
    clearOut();
}


// ************************************************************************* //
