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
#include "pointPatchFields.H"
#include "valuePointPatchFields.H"

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

    if (debug)
    {
        InfoInFunction
            << "Calculating displaced patch"
            << endl;
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
              + this->interpolator().faceToPointInterpolate
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
        const pointPatchVectorField& disp =
            this->mesh_.lookupObject<pointVectorField>
            (
                displacementField_
            ).boundaryField()[this->patch_.index()];

        pointField pdisp(disp.size());
        if (isA<valuePointPatchVectorField>(disp))
        {
            pdisp = dynamicCast<const valuePointPatchVectorField>(disp);
        }
        else
        {
            pdisp = disp.patchInternalField();
        }

        dispPatchPtr_.set
        (
            new standAlonePatch
            (
                this->globalPatch(),
                this->globalPatch().localPoints()
              + this->patchPointToGlobal(pdisp)
            )
        );
    }
    else
    {
        dispPatchPtr_.set(new standAlonePatch(globalPatch()));
    }
    return dispPatchPtr_();
}


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
            dispPatch(),
            samplePatch().dispPatch(),
            intersection::algorithm::visible,
            intersection::direction::contactSphere
        )
    );
}


void Foam::coupledGlobalPolyPatch::clearOut() const
{
    clearInterp();
    globalPolyPatch::clearOut();
}


void Foam::coupledGlobalPolyPatch::clearInterp() const
{
    if
    (
        dispPatchPtr_.valid()
     || patchToSampleInterpPtr_.valid()
    )
    {
        dispPatchPtr_.clear();
        patchToSampleInterpPtr_.clear();
        samplePatch().clearInterp();
    }
    globalPolyPatch::clearInterp();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::coupledGlobalPolyPatch::coupledGlobalPolyPatch
(
    const word& patchName,
    const polyMesh& mesh,
    const word& displacementField
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
    dispPatchPtr_(nullptr),
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
    dispPatchPtr_(nullptr),
    patchToSampleInterpPtr_(nullptr)
{}


Foam::coupledGlobalPolyPatch::coupledGlobalPolyPatch
(
    const polyPatch& pp,
    const word& displacementField
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
    dispPatchPtr_(nullptr),
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
    if (!patchToSampleInterpPtr_.valid())
    {
        calcPatchToPatchInterp();
    }
    return patchToSampleInterpPtr_();
}


void Foam::coupledGlobalPolyPatch::movePoints()
{
    Info<<"cgpp movePoints"<<endl;
    clearOut();
}


void Foam::coupledGlobalPolyPatch::updateMesh()
{
    clearOut();
}


// ************************************************************************* //
