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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(false),
    limitCoeff_(1.0),
    relaxFac_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0)),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check if traction is time-varying
    if (dict.found("tractionSeries"))
    {
        Info<< "    traction is time-varying" << endl;
        tractionSeries_ = Function1<vector>::New("tractionSeries", dict);
        traction_ = tractionSeries_->value(this->db().time().value());
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        Info<< "    pressure is time-varying" << endl;
        pressureSeries_ = Function1<scalar>::New("pressureSeries", dict);
        pressure_ = pressureSeries_->value(this->db().time().value());
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

    if (secondOrder_)
    {
        Info<< "    second order correction" << endl;
    }

    if (limitCoeff_)
    {
        Info<< "    limiter coefficient: " << limitCoeff_ << endl;
    }

    if (relaxFac_ < 1.0)
    {
        Info<< "    relaxation factor: " << relaxFac_ << endl;
    }
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    traction_(mapper(stpvf.traction_)),
    pressure_(mapper(stpvf.pressure_)),
    tractionSeries_(stpvf.tractionSeries_, false),
    pressureSeries_(stpvf.pressureSeries_, false),
    secondOrder_(stpvf.secondOrder_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_, false),
    pressureSeries_(stpvf.pressureSeries_, false),
    secondOrder_(stpvf.secondOrder_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    m(traction_, traction_);
    m(pressure_, pressure_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const solidTractionFvPatchVectorField& dmptf =
        refCast<const solidTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (tractionSeries_.valid())
    {
        traction_ = tractionSeries_->value(this->db().time().value());
    }

    if (pressureSeries_.valid())
    {
        pressure_ = pressureSeries_->value(this->db().time().value());
    }

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    gradient() =
        relaxFac_*solMod.tractionBoundarySnGrad
        (
            traction_, pressure_, patch()
        );
    gradient() += (1.0 - relaxFac_)*gradient();

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Face unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    if (secondOrder_)
    {
        const vectorField dUP(k & gradField.patchInternalField());
        const vectorField nGradUP(n & gradField.patchInternalField());

        Field<vector>::operator=
        (
            patchInternalField()
          + dUP
          + 0.5*(gradient() + nGradUP)/patch().deltaCoeffs()
        );
    }
    else
    {

        Field<vector>::operator=
        (
            patchInternalField()
          + (k & gradField.patchInternalField())
          + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


void solidTractionFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //fixedGradientFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

    if (tractionSeries_.valid())
    {
        writeEntry(os, tractionSeries_());
    }
    else
    {
        writeEntry(os, "traction", traction_);
    }

    if (pressureSeries_.valid())
    {
        writeEntry(os, pressureSeries_());
    }
    else
    {
        writeEntry(os, "pressure", pressure_);
    }
    writeEntry(os, "secondOrder", secondOrder_);
    writeEntry(os, "limitCoeff", limitCoeff_);
    writeEntry(os, "relaxFac", relaxFac_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
