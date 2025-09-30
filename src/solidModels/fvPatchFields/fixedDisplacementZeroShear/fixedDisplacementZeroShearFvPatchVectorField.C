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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedDisplacementZeroShearFvPatchVectorField::
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


Foam::fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& fdzspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(fdzspvf, p, iF, mapper),
    totalDisp_(mapper(fdzspvf.totalDisp_)),
    dispSeries_(fdzspvf.dispSeries_, false),
    forceZeroShearGrad_(fdzspvf.forceZeroShearGrad_)
{}


Foam::fixedDisplacementZeroShearFvPatchVectorField::
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
        DebugInfo
            << "    displacement is time-varying" << endl;
        dispSeries_ =
            Function1<vector>::New
            (
                "displacementSeries",
                this->db().time().userUnits(),
                dimLength,
                dict
            );

        refValue() = dispSeries_->value(this->db().time().value());
    }
    else if (dict.found("value"))
    {
        refValue() = vectorField("value", dict, p.size());
    }
    else
    {
        FatalErrorInFunction
            << "value entry not found for patch " << patch().name()
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


Foam::fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& fdzspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(fdzspvf, iF),
    totalDisp_(fdzspvf.totalDisp_),
    dispSeries_(fdzspvf.dispSeries_, false),
    forceZeroShearGrad_(fdzspvf.forceZeroShearGrad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedDisplacementZeroShearFvPatchVectorField::map
(
    const fvPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    solidDirectionMixedFvPatchVectorField::map(ptf, mapper);

    const fixedDisplacementZeroShearFvPatchVectorField& fdzspvf =
        refCast<const fixedDisplacementZeroShearFvPatchVectorField>(ptf);

    mapper(totalDisp_, fdzspvf.totalDisp_);
}


void Foam::fixedDisplacementZeroShearFvPatchVectorField::reset
(
    const fvPatchField<vector>& ptf
)
{
    solidDirectionMixedFvPatchVectorField::reset(ptf);

    const fixedDisplacementZeroShearFvPatchVectorField& fdzspvf =
        refCast<const fixedDisplacementZeroShearFvPatchVectorField>(ptf);

    totalDisp_.reset(fdzspvf.totalDisp_);
}


void Foam::fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp(totalDisp_);

    if (dispSeries_.valid())
    {
        disp = dispSeries_->value(this->db().time().value());
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


void Foam::fixedDisplacementZeroShearFvPatchVectorField::write
(
    Ostream& os
) const
{
    if (dispSeries_.valid())
    {
        writeEntry(os, "displacementSeries", dispSeries_());
    }
    writeEntry(os, "forceZeroShearGrad", forceZeroShearGrad_);

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fixedDisplacementZeroShearFvPatchVectorField
    );
}

// ************************************************************************* //
