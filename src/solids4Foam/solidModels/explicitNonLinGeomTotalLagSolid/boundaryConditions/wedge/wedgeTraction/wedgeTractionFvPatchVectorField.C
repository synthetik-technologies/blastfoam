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

#include "wedgeTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wedgeTractionFvPatchVectorField::wedgeTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    wedgeFvPatchVectorField(p, iF)
{}


wedgeTractionFvPatchVectorField::wedgeTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    wedgeFvPatchVectorField(p, iF)
{}


wedgeTractionFvPatchVectorField::wedgeTractionFvPatchVectorField
(
    const wedgeTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wedgeFvPatchVectorField(ptf, p, iF, mapper)
{}


wedgeTractionFvPatchVectorField::wedgeTractionFvPatchVectorField
(
    const wedgeTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    wedgeFvPatchVectorField(rifvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wedgeTractionFvPatchVectorField::evaluate
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
            "grad(rhoU)"
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

    const fvsPatchField<vector>& n =
        patch().lookupPatchField<surfaceVectorField, vector>("n");

    const Field<tensor> gradPx
    (
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(P.x)"
        ).patchInternalField()
    );
    const Field<tensor> gradPy
    (
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(P.y)"
        ).patchInternalField()
    );
    const Field<tensor> gradPz
    (
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(P.z)"
        ).patchInternalField()
    );

    Field<tensor> gradT(this->size());
    forAll(gradT, i)
    {
        gradT[i] = tensor
        (
            (gradPx[i] & n[i]),
            (gradPy[i] & n[i]),
            (gradPz[i] & n[i])
        );
    }
    Field<vector> pif(this->patchInternalField());
    pif += (k & gradT);

    Field<vector>::operator=
    (
        transform(wedgePatch.faceT(), pif)
    );
}


void wedgeTractionFvPatchVectorField::write(Ostream& os) const
{
    wedgeFvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    wedgeTractionFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
