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

#include "thermalConvectionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalConvectionFvPatchScalarField::thermalConvectionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    DTName_("undefined"),
    alpha_(p.size(), 0),
    Tinf_("Tinf", dimTemperature, 293)
{}


Foam::thermalConvectionFvPatchScalarField::thermalConvectionFvPatchScalarField
(
    const thermalConvectionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    DTName_(ptf.DTName_),
    alpha_(mapper(ptf.alpha_)),
    Tinf_(ptf.Tinf_)
{}


Foam::thermalConvectionFvPatchScalarField::thermalConvectionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    DTName_(dict.lookupOrDefault<word>("thermalConductivityName", "k")),
    alpha_("alpha", dict, p.size()),
    Tinf_("Tinf", dimTemperature, readScalar(dict.lookup("Tinf")))
{
    Info<< patch().name() << ": thermalConvection" << endl;

    fvPatchField<scalar>::operator=(patchInternalField());
}


Foam::thermalConvectionFvPatchScalarField::thermalConvectionFvPatchScalarField
(
    const thermalConvectionFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wbppsf, iF),
    DTName_(wbppsf.DTName_),
    alpha_(wbppsf.alpha_),
    Tinf_(wbppsf.Tinf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::thermalConvectionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(alpha_, alpha_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::thermalConvectionFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const thermalConvectionFvPatchScalarField& tiptf =
        refCast<const thermalConvectionFvPatchScalarField>(ptf);

    //alpha_.resize(tiptf.alpha_.size());
    alpha_.rmap(tiptf.alpha_, addr);
}


void Foam::thermalConvectionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);
    scalarField TP(this->patchInternalField());

    // Lookup thermal diffusivity i.e. conductivity
    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

    Field<scalar>::operator=
    (
        (alpha_*Tinf_.value() + DT*TP/delta)
       /(DT/delta + alpha_ + SMALL)
    );

    fvPatchField<scalar>::evaluate();
}

Foam::tmp<Foam::Field<Foam::scalar> >
Foam::thermalConvectionFvPatchScalarField::snGrad() const
{
    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);
    scalarField TP(this->patchInternalField());

    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

//     scalarField gradient =
//          DT*alpha_*(Tinf_.value() - TP)
//         /(DT + alpha_*delta + SMALL);

//     scalarField gradient =
//          alpha_*(Tinf_.value() - TP)
//         /(DT + alpha_*delta + SMALL);

    return tmp<Field<scalar> >
    (
         alpha_*(Tinf_.value() - TP)
        /(DT + alpha_*delta + SMALL)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::thermalConvectionFvPatchScalarField::gradientInternalCoeffs() const
{
    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);

    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

//     return tmp<Field<scalar> >
//     (
//         -pTraits<scalar>::one*DT*alpha_/(DT + alpha_*delta + SMALL)
//     );

    return tmp<Field<scalar> >
    (
        -pTraits<scalar>::one*alpha_/(DT + alpha_*delta + SMALL)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::thermalConvectionFvPatchScalarField::gradientBoundaryCoeffs() const
{
    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);

    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

    return alpha_*Tinf_.value()/(DT + alpha_*delta + SMALL);
}


void Foam::thermalConvectionFvPatchScalarField::write(Ostream& os) const
{
    writeEntry(os, "thermalConductivityName", DTName_);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "Tinf", Tinf_);

    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        thermalConvectionFvPatchScalarField
    );
}

// ************************************************************************* //
