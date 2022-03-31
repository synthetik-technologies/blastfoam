/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solidSymmetryFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(p, iF),
    secondOrder_(false)
{}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    symmetryFvPatchField<vector>(ptf, p, iF, mapper),
    secondOrder_(ptf.secondOrder_)
{
    if (!isType<symmetryFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "solidSymmetryFvPatchVectorField::"
            "solidSymmetryFvPatchVectorField\n"
            "(\n"
            "    const solidSymmetryFvPatchVectorField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    symmetryFvPatchField<vector>(p, iF, dict),
    secondOrder_(false)
{
    Info << "Symmetry boundary condition with non-orthogonal correction"
        << endl;

    if (dict.found("secondOrder"))
    {
        secondOrder_ = Switch(dict.lookup("secondOrder"));
        Info<< "Second order correction: " << secondOrder_ << endl;
    }

    if (!isType<symmetryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "solidSymmetryFvPatchVectorField::"
            "solidSymmetryFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<vector>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(ptf, iF),
    secondOrder_(ptf.secondOrder_)
{}


// return gradient at boundary
tmp<Field<vector> > solidSymmetryFvPatchVectorField::snGrad() const
{
    // Unit normals
    const vectorField nHat(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(nHat)) & delta);

    // Lookup the gradient of displacement field
    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Calculate the corrected patch internal field
    const vectorField DP
    (
        patchInternalField()
      + (k & gradD.patchInternalField())
    );

    if (secondOrder_)
    {
        // Normal component of patch internal gradient
        const vectorField nGradDP(nHat & gradD.patchInternalField());

        return
          2*(
                transform(I - 2.0*sqr(nHat), DP) - DP
            )*(patch().deltaCoeffs()/2.0)
          - transform(sqr(nHat), nGradDP);
    }
    else
    {
        return
        (
            transform(I - 2.0*sqr(nHat), DP)
          - DP
        )*(patch().deltaCoeffs()/2.0);
    }
}


// Evaluate the field on the patch
void solidSymmetryFvPatchVectorField::
evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Unit normals
    const vectorField nHat(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(nHat)) & delta);

    // Lookup the gradient of displacement field
    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Calculate the corrected patch internal field
    const vectorField DP
    (
        patchInternalField()
      + (k & gradD.patchInternalField())
    );

    if (secondOrder_)
    {
        vectorField nGradDP(nHat&gradD.patchInternalField());

        Field<vector>::operator=
        (
            transform
            (
                I - sqr(nHat),
                DP + 0.5*nGradDP/patch().deltaCoeffs()
            )
        );
    }
    else
    {
        Field<vector>::operator=
        (
            (
                DP
              + transform(I - 2.0*sqr(nHat), DP)
            )/2.0
        );
    }

    transformFvPatchField<vector>::evaluate();
}


// Write
void solidSymmetryFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "secondOrder", secondOrder_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidSymmetryFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
