/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technologies: Mapping of patches
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "mappedMovingPatchBase.H"
#include "volFields.H"
#include "pointFields.H"
#include "globalPolyBoundaryMesh.H"
#include "coupledGlobalPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedMovingPatchBase, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mappedMovingPatchBase::setOffsets() const
{
    const polyMesh& mesh = this->patch_.boundaryMesh().mesh();

    vectorField& offsetsRef = const_cast<vectorField&>(this->offsets_);
    offsetsRef.resize(this->patch_.size());
    offsetsRef = vector::zero;
    if (!returnReduce(this->patch_.size(), maxOp<label>()))
    {
        return;
    }
    if (mesh.foundObject<volVectorField>(displacementField_))
    {
        offsetsRef +=
            mesh.lookupObject<volVectorField>
            (
                displacementField_
            ).boundaryField()[this->patch_.index()];
    }
    else if (mesh.foundObject<pointVectorField>(displacementField_))
    {
        offsetsRef +=
            mesh.lookupObject<pointVectorField>
            (
                displacementField_
            ).boundaryField()[this->patch_.index()].patchInternalField();
    }

    if
    (
        !this->patch_.boundaryMesh().mesh().time().foundObject<polyMesh>
        (
            sampleRegion()
        )
    )
    {
        return;
    }

    const polyPatch& sPatch(this->samplePolyPatch());
    if (!returnReduce(sPatch.size(), maxOp<label>()))
    {
        return;
    }

    const polyMesh& sMesh(this->sampleMesh());
    if (sMesh.foundObject<volVectorField>(sampleDisplacementField_))
    {
        const coupledGlobalPolyPatch& cgpp =
            dynamicCast<const coupledGlobalPolyPatch>
            (
                globalPolyBoundaryMesh::New(sMesh)[sPatch]
            );
        offsetsRef -= cgpp.faceInterpolate
        (
            sMesh.lookupObject<volVectorField>
            (
                sampleDisplacementField_
            ).boundaryField()[sPatch.index()]
        );
    }
    else if (sMesh.foundObject<pointVectorField>(sampleDisplacementField_))
    {
        const coupledGlobalPolyPatch& cgpp =
            dynamicCast<const coupledGlobalPolyPatch>
            (
                globalPolyBoundaryMesh::New(sMesh)[sPatch]
            );
        offsetsRef -= cgpp.pointToFaceInterpolate
        (
            sMesh.lookupObject<pointVectorField>
            (
                sampleDisplacementField_
            ).boundaryField()[sPatch.index()].patchInternalField()
        );
    }
}


void Foam::mappedMovingPatchBase::calcMapping() const
{
    setOffsets();
    mappedPatchBase::calcMapping();
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp
)
:
    mappedPatchBase(pp),
    displacementField_(word::null),
    sampleDisplacementField_(word::null)
{
    this->offsetMode_ = mappedPatchBase::NONUNIFORM;
}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const mappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offsets
)
:
    mappedPatchBase
    (
        pp,
        sampleRegion,
        mode,
        samplePatch,
        offsets
    ),
    displacementField_(word::null),
    sampleDisplacementField_(word::null)
{
    this->offsetMode_ = mappedPatchBase::NONUNIFORM;
}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    mappedPatchBase
    (
        pp,
        dict.lookupOrDefault<word>("sampleRegion", ""),
        this->sampleModeNames_.read(dict.lookup("sampleMode")),
        dict.lookupOrDefault<word>("samplePatch", ""),
        vectorField(pp.size(), Zero)
    ),
    displacementField_
    (
        dict.lookupOrDefault("displacementField", word::null)
    ),
    sampleDisplacementField_
    (
        dict.lookupOrDefault("sampleDisplacementField", word::null)
    )
{
    this->offsetMode_ = mappedPatchBase::NONUNIFORM;
}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const mappedMovingPatchBase& mmpb
)
:
    mappedPatchBase(pp, mmpb),
    displacementField_(mmpb.displacementField_),
    sampleDisplacementField_(mmpb.sampleDisplacementField_)
{
    this->offsetMode_ = mappedPatchBase::NONUNIFORM;
}


Foam::mappedMovingPatchBase::mappedMovingPatchBase
(
    const polyPatch& pp,
    const mappedMovingPatchBase& mmpb,
    const labelUList& mapAddressing
)
:
    mappedPatchBase(pp, mmpb, mapAddressing),
    displacementField_(mmpb.displacementField_),
    sampleDisplacementField_(mmpb.sampleDisplacementField_)
{
    this->offsetMode_ = mappedPatchBase::NONUNIFORM;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedMovingPatchBase::~mappedMovingPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedMovingPatchBase::write(Ostream& os) const
{
    mappedPatchBase::write(os);
    writeEntryIfDifferent
    (
        os,
        "displacementField",
        word::null,
        displacementField_
    );
    writeEntryIfDifferent
    (
        os,
        "sampleDisplacementField",
        word::null,
        sampleDisplacementField_
    );
}
// ************************************************************************* //
