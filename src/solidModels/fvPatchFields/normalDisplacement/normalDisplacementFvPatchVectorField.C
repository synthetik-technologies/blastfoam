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

#include "normalDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalDisplacementFvPatchVectorField::
normalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    normalDisp_(p.size(), 0.0),
    dispSeries_()
{}


Foam::normalDisplacementFvPatchVectorField::
normalDisplacementFvPatchVectorField
(
    const normalDisplacementFvPatchVectorField& ndpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ndpvf, p, iF, mapper),
    normalDisp_(mapper(ndpvf.normalDisp_)),
    dispSeries_(ndpvf.dispSeries_, false)
{}


Foam::normalDisplacementFvPatchVectorField::
normalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    normalDisp_(p.size(), 0.0),
    dispSeries_()
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        DebugInfo
            << "    normal displacement is time-varying" << endl;
        dispSeries_ = Function1<scalar>::New
        (
            "displacementSeries",
            this->db().time().userUnits(),
            dimLength,
            dict
        );

        fvPatchField<vector>::operator==
        (
            patch().nf()*dispSeries_->value(this->db().time().value())
        );
    }
    else
    {
        normalDisp_ = scalarField("normalDisplacement", dict, p.size());

        fvPatchField<vector>::operator==
        (
            patch().nf()*normalDisp_
        );
    }
}


Foam::normalDisplacementFvPatchVectorField::
normalDisplacementFvPatchVectorField
(
    const normalDisplacementFvPatchVectorField& ndpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ndpvf, iF),
    normalDisp_(ndpvf.normalDisp_),
    dispSeries_(ndpvf.dispSeries_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::normalDisplacementFvPatchVectorField::map
(
    const fvPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    const normalDisplacementFvPatchVectorField& ndpvf =
        refCast<const normalDisplacementFvPatchVectorField>(ptf);

    mapper(normalDisp_, ndpvf.normalDisp_);
}


void Foam::normalDisplacementFvPatchVectorField::reset
(
    const fvPatchField<vector>& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    const normalDisplacementFvPatchVectorField& ndpvf =
        refCast<const normalDisplacementFvPatchVectorField>(ptf);

    normalDisp_.reset(ndpvf.normalDisp_);
}


void Foam::normalDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField nDisp = normalDisp_;

    if (dispSeries_.valid())
    {
        nDisp = dispSeries_->value(this->db().time().value());
    }

    vectorField disp(nDisp*patch().nf());

    if (internalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    fvPatchField<vector>::operator==(disp);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::normalDisplacementFvPatchVectorField::snGrad() const
{
    //- fixedValue snGrad with no correction
    //  return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Unit normal vectors
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vectors
    const vectorField k((I - sqr(n)) & delta);

    return
    (
        *this
        - (patchInternalField() + (k & gradField.patchInternalField()))
    )*patch().deltaCoeffs();
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::normalDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Unit normal vectors
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vectors
    const vectorField k((I - sqr(n)) & delta);

    return
    (
        this->patch().deltaCoeffs()
       *(*this - (k & gradField.patchInternalField()))
    );
}


void Foam::normalDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.valid())
    {
        writeEntry(os, "displacementSeries", dispSeries_());
    }
    else
    {
        writeEntry(os, "normalDisp", normalDisp_);
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        normalDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
