/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "extrapolatedPressureValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedPressureValueFvPatchScalarField::extrapolatedPressureValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    relaxFac_(1.0)
{}


extrapolatedPressureValueFvPatchScalarField::extrapolatedPressureValueFvPatchScalarField
(
    const extrapolatedPressureValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    relaxFac_(1.0)
{}


extrapolatedPressureValueFvPatchScalarField::extrapolatedPressureValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    relaxFac_(readScalar(dict.lookup("relaxFactor")))
{}


extrapolatedPressureValueFvPatchScalarField::extrapolatedPressureValueFvPatchScalarField
(
    const extrapolatedPressureValueFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    relaxFac_(pivpvf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedPressureValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Name of the gradient field
    const word gradFieldName("grad(" + dimensionedInternalField().name() + ")");

    // If available, use the gradient for non-orthogonal correction
    if (db().foundObject<volVectorField>(gradFieldName))
    {
        // Lookup gradient of the field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>(gradFieldName);

        // Unit normals
        const vectorField n = patch().nf();

        // Delta vectors
        const vectorField delta = patch().delta();

        // Linearly extrapolate from the internal field and apply
        // under-relaxation
        fvPatchField<scalar>::operator==
        (
            (1.0 - relaxFac_)*(*this)
          + relaxFac_
           *(
               patchInternalField() + (delta & gradField.patchInternalField())
            )
        );
    }
    else
    {
        // Zero gradient
        fvPatchField<scalar>::operator==
        (
            (1.0 - relaxFac_)*(*this) + relaxFac_*patchInternalField()
        );
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<Foam::Field<scalar> > extrapolatedPressureValueFvPatchScalarField::snGrad() const
{
    // Name of the gradient field
    const word gradFieldName("grad(" + dimensionedInternalField().name() + ")");

    // If available, use the gradient for non-orthogonal correction
    if (db().foundObject<volVectorField>(gradFieldName))
    {
        // Lookup gradient of the field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>(gradFieldName);

        // Unit normals
        const vectorField n = patch().nf();

        // Delta vectors
        const vectorField delta = patch().delta();

        // Correction vectors
        const vectorField k = (I - sqr(n)) & delta;

        // Surface normal gradient using non-orthogonal correction
        return
        (
            *this - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        // Surface normal gradient without non-orthogonal correction
        return
        (
            *this - patchInternalField()
        )*patch().deltaCoeffs();
    }
}

tmp<Field<scalar> > extrapolatedPressureValueFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    // Name of the gradient field
    const word gradFieldName("grad(" + dimensionedInternalField().name() + ")");

    // If available, use the gradient for non-orthogonal correction
    if (db().foundObject<volVectorField>(gradFieldName))
    {
        // Lookup gradoent of the field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>(gradFieldName);

        // Unit normals
        const vectorField n = patch().nf();

        // Delta vectors
        const vectorField delta = patch().delta();

        // Correction vectors
        const vectorField k = (I - sqr(n)) & delta;

        // Coefficients with non-orthogonal correction
        return patch().deltaCoeffs()*(*this - (k & gradField.patchInternalField()));
    }
    else
    {
        // Coefficients without non-orthogonal correction
        return patch().deltaCoeffs()*(*this);
    }
}

void extrapolatedPressureValueFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("relaxFactor")
        << relaxFac_ << token::END_STATEMENT << nl;

    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    extrapolatedPressureValueFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
