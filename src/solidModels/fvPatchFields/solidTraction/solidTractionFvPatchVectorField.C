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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    tractionBase(p),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(false),
    relaxFac_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


Foam::solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool mustRead
)
:
    fixedGradientFvPatchVectorField(p, iF),
    tractionBase(p, dict, false),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    relaxFac_()
{
    DebugInfo
        << "Creating " << type() << " boundary condition" << endl;

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
        DebugInfo<< "    traction is time-varying" << endl;
        tractionSeries_ = Function1<vector>::New
        (
            "tractionSeries",
            this->db().time().userUnits(),
            dimPressure,
            dict
        );
        this->traction() =
            tractionSeries_->value(this->db().time().value());
    }
    else if (mustRead)
    {
        this->traction() = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        DebugInfo<< "    pressure is time-varying" << endl;
        pressureSeries_ = Function1<scalar>::New
        (
            "pressureSeries",
            this->db().time().userUnits(),
            dimPressure,
            dict
        );
        this->pressure() =
            pressureSeries_->value(this->db().time().value());
    }
    else if (mustRead)
    {
        this->pressure() = scalarField("pressure", dict, p.size());
    }

    if (dict.found("relaxationFactor"))
    {
        DebugInfo<< "    Using relaxationFactor" << endl;
        relaxFac_ = Function1<scalar>::New
        (
            "relaxationFactor",
            this->db().time().userUnits(),
            dimless,
            dict
        );
    }

    if (secondOrder_)
    {
        DebugInfo<< "    second order correction" << endl;
    }
}


Foam::solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    tractionBase(stpvf, p, mapper),
    tractionSeries_(stpvf.tractionSeries_, false),
    pressureSeries_(stpvf.pressureSeries_, false),
    secondOrder_(stpvf.secondOrder_),
    relaxFac_(stpvf.relaxFac_, false)
{}


Foam::solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    tractionBase(stpvf, this->patch()),
    tractionSeries_(stpvf.tractionSeries_, false),
    pressureSeries_(stpvf.pressureSeries_, false),
    secondOrder_(stpvf.secondOrder_),
    relaxFac_(stpvf.relaxFac_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidTractionFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    fixedGradientFvPatchVectorField::map(ptf, mapper);
    tractionBase::map(ptf, mapper);
}


void Foam::solidTractionFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedGradientFvPatchVectorField::reset(ptf);
    tractionBase::reset(ptf);
}


bool Foam::solidTractionFvPatchVectorField::updateFields()
{
    bool updateTraction = false;
    if (tractionSeries_.valid())
    {
        this->traction() =
            tractionSeries_->value(this->db().time().value());
        updateTraction = true;
    }

    if (pressureSeries_.valid())
    {
        this->pressure() =
            pressureSeries_->value(this->db().time().value());
    }

    updateForce();

    return updateTraction;
}


void Foam::solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateFields();

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    if (relaxFac_.valid() && canRelax)
    {
        scalar relaxFac = relaxFac_->value(this->db().time().value());
        gradient() =
            relaxFac*solMod.tractionBoundarySnGrad
            (
                this->traction(), this->pressure(), patch()
            )
          + (1.0 - relaxFac)*gradient();
    }
    else
    {
        gradient() =
            solMod.tractionBoundarySnGrad
            (
                this->traction(), this->pressure(), patch()
            );
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::solidTractionFvPatchVectorField::evaluate
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
    const vectorField n(this->patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((tensor::I - sqr(n)) & delta);

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


void Foam::solidTractionFvPatchVectorField::write(Ostream& os) const
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
        writeEntry(os, "traction", this->traction());
    }

    if (pressureSeries_.valid())
    {
        writeEntry(os, pressureSeries_());
    }
    else
    {
        writeEntry(os, "pressure", this->pressure());
    }

    if (relaxFac_.valid())
    {
        writeEntry(os, relaxFac_());
    }

    writeEntry(os, "secondOrder", secondOrder_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        solidTractionFvPatchVectorField
    );
}


// ************************************************************************* //
