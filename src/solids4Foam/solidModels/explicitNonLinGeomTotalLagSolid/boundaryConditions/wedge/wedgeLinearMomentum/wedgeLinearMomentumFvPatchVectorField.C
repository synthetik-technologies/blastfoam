/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "wedgeLinearMomentumFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wedgeLinearMomentumFvPatchVectorField::
wedgeLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    wedgeFvPatchVectorField(p, iF)
{}


wedgeLinearMomentumFvPatchVectorField::
wedgeLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    wedgeFvPatchVectorField(p, iF, dict)
{}


wedgeLinearMomentumFvPatchVectorField::
wedgeLinearMomentumFvPatchVectorField
(
    const wedgeLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wedgeFvPatchVectorField(ptf, p, iF, mapper)
{}


wedgeLinearMomentumFvPatchVectorField::
wedgeLinearMomentumFvPatchVectorField
(
    const wedgeLinearMomentumFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    wedgeFvPatchVectorField(rifvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wedgeLinearMomentumFvPatchVectorField::evaluate
(
    const Pstream::commsTypes comms
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if
    (
        !patch().boundaryMesh().mesh().foundObject<volTensorField>
        (
            "grad(" + this->internalField().name() + ")"
        )
    )
    {
        wedgeFvPatchVectorField::evaluate(comms);
        return;
    }
    const wedgeFvPatch& wedgePatch =
        refCast<const wedgeFvPatch>(this->patch());

    // Rotate patchC field back to centre plane to find transformed cell centres
    const vectorField patchC(patch().patch().faceCentres());
    vectorField transC(wedgePatch.faceT().T() & patchC);

    // Calculate correction vector which connects actual cell centre to the
    // transformed cell centre
    const vectorField k(transC - patch().Cn());

    const fvPatchField<tensor>& gradLM =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->internalField().name() + ")"
        );

    Field<vector> pif(this->patchInternalField());
    pif += (k & gradLM.patchInternalField());

    Field<vector>::operator=
    (
        transform(wedgePatch.faceT(), pif)
    );
}


void wedgeLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    wedgeFvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    wedgeLinearMomentumFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
