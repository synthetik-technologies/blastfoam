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

#include "holeVelocityCyclicFvPatchVectorField.H"
#include "cyclicFvPatchFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    holePatchField(p),
    cyclicPatchName_(),
    cyclicPatchLabel_(-1),
    initWallSf_(0),
    initCyclicSf_(0),
    nbrCyclicSf_(0),
    curTimeIndex_(-1)
{}


Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    holePatchField(p, dict),
    cyclicPatchName_(dict.lookup("cyclicPatch")),
    cyclicPatchLabel_(p.patch().boundaryMesh().findPatchID(cyclicPatchName_)),
    initWallSf_(0),
    initCyclicSf_(0),
    nbrCyclicSf_(0),
    curTimeIndex_(-1)
{
    if (p.size() > 0)
    {
        initWallSf_ = p.Sf();
        initCyclicSf_ = p.boundaryMesh()[cyclicPatchLabel_].Sf();
        nbrCyclicSf_ =  refCast<const cyclicFvPatch>
        (
            p.boundaryMesh()[cyclicPatchLabel_]
        ).neighbFvPatch().Sf();
    }
}


Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const holeVelocityCyclicFvPatchVectorField& hpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(hpf, p, iF, mapper),
    holePatchField(hpf, p, mapper),
    cyclicPatchName_(hpf.cyclicPatchName_),
    cyclicPatchLabel_(hpf.cyclicPatchLabel_),
    initWallSf_(hpf.initWallSf_),
    initCyclicSf_(hpf.initCyclicSf_),
    nbrCyclicSf_(hpf.nbrCyclicSf_),
    curTimeIndex_(-1)
{}


Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const holeVelocityCyclicFvPatchVectorField& hpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(hpf, iF),
    holePatchField(hpf),
    cyclicPatchName_(hpf.cyclicPatchName_),
    cyclicPatchLabel_(hpf.cyclicPatchLabel_),
    initWallSf_(hpf.initWallSf_),
    initCyclicSf_(hpf.initCyclicSf_),
    nbrCyclicSf_(hpf.nbrCyclicSf_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::holeVelocityCyclicFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    holePatchField::autoMap(m);

    if (patch().size() > 0)
    {
        const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
        initWallSf_ = patch().patchSlice(areas);
        initCyclicSf_ = patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ].patchSlice(areas);
        nbrCyclicSf_ = refCast<const cyclicFvPatch>
        (
            patch().boundaryMesh()
            [
                cyclicPatchLabel_
            ]
        ).neighbFvPatch().patch().patchSlice(areas);
    }

}


void Foam::holeVelocityCyclicFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const holeVelocityCyclicFvPatchVectorField& hpf =
        refCast<const holeVelocityCyclicFvPatchVectorField>(ptf);
    holePatchField::rmap(hpf, addr);

    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::holeVelocityCyclicFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    holePatchField::updateCoeffs();

    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const fvPatch& cyclicPatch =
            patch().boundaryMesh()[cyclicPatchLabel_];
        const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
        (
            cyclicPatch
        ).neighbFvPatch();

        vectorField::subField Sfw = patch().patch().faceAreas();
        vectorField newSfw((1.0 - openFraction_)*initWallSf_);
        forAll(Sfw, facei)
        {
            Sfw[facei] = newSfw[facei];
        }
        const_cast<scalarField&>(patch().magSf()) = mag(patch().Sf());

        // Update owner side of cyclic
        const_cast<vectorField&>(cyclicPatch.Sf()) =
            openFraction_*initCyclicSf_;

        const_cast<scalarField&>(cyclicPatch.magSf()) =
            mag(cyclicPatch.Sf());

        // Update neighbour side of cyclic
        const_cast<vectorField&>(nbrPatch.Sf()) =
            openFraction_*nbrCyclicSf_;

        const_cast<scalarField&>(nbrPatch.magSf()) =
            mag(nbrPatch.Sf());

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::holeVelocityCyclicFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    writeEntry(os, "cyclicPatch", cyclicPatchName_);
    writeEntry(os, "L", L_);
    writeEntry(os, "LName", LName_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        holeVelocityCyclicFvPatchVectorField
    );
}

// ************************************************************************* //
