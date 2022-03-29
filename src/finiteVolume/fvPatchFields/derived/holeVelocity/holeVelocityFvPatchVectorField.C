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
#include "noSlipFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    cyclicFvPatchVectorField(p, iF),
    L_(0.0),
    LName_(IOobject::groupName("d", iF.group())),
    mask_(p.size(), 0.0),
    wallPatchField_
    (
        new noSlipFvPatchVectorField
        (
            p,
            iF
        )
    )
{}


Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchVectorField(p, iF, dict),
    L_(dict.lookup<scalar>("L")),
    LName_(dict.lookup("LName")),
    mask_(p.size(), 0.0),
    wallPatchField_()
{
    // Create a new patch dictionary and replace the type with the intactType
    dictionary wallDict(dict.parent(), dict.subDict("wallPatch"));
    if (!wallDict.found("patchType"))
    {
        wallDict.add("patchType", typeName);
    }
    wallPatchField_ =
        fvPatchField<vector>::New
        (
            p,
            iF,
            wallDict
        );
}


Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const holeVelocityCyclicFvPatchVectorField& bpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvPatchVectorField(bpf, p, iF, mapper),
    L_(bpf.L_),
    LName_(bpf.LName_),
    mask_(mapper(bpf.mask_)),
    wallPatchField_
    (
        fvPatchField<vector>::New
        (
            bpf.wallPatchField_(),
            p,
            iF,
            mapper
        )
    )
{}


Foam::holeVelocityCyclicFvPatchVectorField::holeVelocityCyclicFvPatchVectorField
(
    const holeVelocityCyclicFvPatchVectorField& bpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    cyclicFvPatchVectorField(bpf, iF),
    L_(bpf.L_),
    LName_(bpf.LName_),
    mask_(bpf.mask_),
    wallPatchField_(bpf.wallPatchField_->clone(iF))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::holeVelocityCyclicFvPatchVectorField::patchNeighbourField() const
{
    return
        mask_*cyclicFvPatchVectorField::patchNeighbourField()
      + (1.0 - mask_)*wallPatchField_();

}


void Foam::holeVelocityCyclicFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    cyclicFvPatchVectorField::autoMap(m);
    wallPatchField_->autoMap(m);
     m(mask_, mask_);
}


void Foam::holeVelocityCyclicFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    cyclicFvPatchVectorField::rmap(ptf, addr);

    const holeVelocityCyclicFvPatchVectorField& bpf =
        refCast<const holeVelocityCyclicFvPatchVectorField>(ptf);
    wallPatchField_->rmap(bpf.wallPatchField_(), addr);
    mask_.rmap(bpf.mask_, addr);
}


void Foam::holeVelocityCyclicFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatchField<scalar>& l =
            patch().lookupPatchField<volScalarField, scalar>(LName_);

    mask_ = pos0(L_ - l);
    wallPatchField_->updateCoeffs();
    cyclicFvPatchVectorField::updateCoeffs();
}


void Foam::holeVelocityCyclicFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    wallPatchField_->evaluate(commsType);
    cyclicFvPatchVectorField::evaluate(commsType);

    vectorField::operator=
    (
        mask_*cyclicFvPatchVectorField::patchNeighbourField()
      + (1.0 - mask_)*wallPatchField_()
    );
}


void Foam::holeVelocityCyclicFvPatchVectorField::write(Ostream& os) const
{
    cyclicFvPatchVectorField::write(os);
    writeKeyword(os, "wallPatch")
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;
    wallPatchField_->write(os);
    os << decrIndent << indent << token::END_BLOCK << endl;
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
