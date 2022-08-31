/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_(),
    forceZeroShearGrad_(false)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    totalDisp_(mapper(ptf.totalDisp_)),
    dispSeries_(ptf.dispSeries_, false),
    forceZeroShearGrad_(ptf.forceZeroShearGrad_)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    totalDisp_("value", dict, p.size()),
    dispSeries_(),
    forceZeroShearGrad_
    (
        dict.lookupOrDefault<Switch>("forceZeroShearGrad", false)
    )
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            Function1<vector>::New("displacementSeries", dict);

        refValue() = dispSeries_->value(this->db().time().timeOutputValue());
    }
    else if (dict.found("value"))
    {
        refValue() = vectorField("value", dict, p.size());
    }
    else
    {
        FatalErrorIn
        (
            "fixedDisplacementZeroShearFvPatchVectorField::"
            "fixedDisplacementZeroShearFvPatchVectorField"
        )   << "value entry not found for patch " << patch().name()
            << abort(FatalError);
    }

    this->refGrad() = vector::zero;

    this->valueFraction() = sqr(patch().nf());

    Field<vector> normalValue(transform(valueFraction(), refValue()));

    Field<vector> gradValue
    (
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs()
    );

    Field<vector> transformGradValue
    (
        transform(I - valueFraction(), gradValue)
    );

    Field<vector>::operator=(normalValue + transformGradValue);
}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    totalDisp_(ptf.totalDisp_),
    dispSeries_(ptf.dispSeries_),
    forceZeroShearGrad_(ptf.forceZeroShearGrad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);

    m(totalDisp_, totalDisp_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const fixedDisplacementZeroShearFvPatchVectorField& dmptf =
        refCast<const fixedDisplacementZeroShearFvPatchVectorField>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp(totalDisp_);

    if (dispSeries_.valid())
    {
        disp = dispSeries_->value(this->db().time().timeOutputValue());
    }

    if (internalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    // Set displacement
    refValue() = disp;

    // Set gradient to zero to force zero shear traction
    if (forceZeroShearGrad_)
    {
        refGrad() = vector::zero;
    }
    else
    {
        // Calculate the shear gradient such that the shear traction is zero

        // Lookup the solidModel object
        const solidModel& solMod =
            lookupSolidModel(patch().boundaryMesh().mesh());

        // Set gradient to force zero shear traction
        refGrad() =
            solMod.tractionBoundarySnGrad
            (
                vectorField(patch().size(), vector::zero),
                scalarField(patch().size(), 0.0),
                patch()
            );
    }

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.valid())
    {
        writeEntry(os, "displacementSeries", dispSeries_());
    }
    writeEntry(os, "forceZeroShearGrad", forceZeroShearGrad_);

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
