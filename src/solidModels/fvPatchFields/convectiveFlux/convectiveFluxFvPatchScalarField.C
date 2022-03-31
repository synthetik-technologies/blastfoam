/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "convectiveFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convectiveFluxFvPatchScalarField::convectiveFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    alpha_(p.size(), 0),
//     lambda_(p.size(), 0),
    Tinf_("Tinf", dimTemperature, 293)
    // UName_("U"),
    // phiName_("phi"),
    // rhoName_("rho"),
    // adjoint_(false)
{}


Foam::convectiveFluxFvPatchScalarField::convectiveFluxFvPatchScalarField
(
    const convectiveFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    alpha_(mapper(ptf.alpha_)),
//     lambda_(ptf.lambda_),
    Tinf_(ptf.Tinf_)
    // UName_(ptf.UName_),
    // phiName_(ptf.phiName_),
    // rhoName_(ptf.rhoName_),
    // adjoint_(ptf.adjoint_)
{}


Foam::convectiveFluxFvPatchScalarField::convectiveFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    alpha_("alpha", dict, p.size()),
//     lambda_(p.size(), 0),
    Tinf_(dict.lookup("Tinf"))
    // UName_(dict.lookupOrDefault<word>("U", "U")),
    // phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    // rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    // adjoint_(dict.lookup("adjoint"))
{
    fvPatchField<scalar>::operator=(patchInternalField());

    // const & lambda =
    //     db().lookupObject<volScalarField>("lambda");

//     const fvPatchField<scalar>& lambda =
//         patch().lookupPatchField<volScalarField, scalar>("lambdaEff");

//     lambda_ = lambda;

    // const dictionary& transportProperties =
    //     db().lookupObject<IOdictionary>("transportProperties");

    // lambda_ = dimensionedScalar(transportProperties.lookup("lambda")).value();
}


Foam::convectiveFluxFvPatchScalarField::convectiveFluxFvPatchScalarField
(
    const convectiveFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wbppsf, iF),
    alpha_(wbppsf.alpha_),
//     lambda_(wbppsf.lambda_),
    Tinf_(wbppsf.Tinf_)
    // UName_(wbppsf.UName_),
    // phiName_(wbppsf.phiName_),
    // rhoName_(wbppsf.rhoName_),
    // adjoint_(wbppsf.adjoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::convectiveFluxFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);
    scalarField TP(this->patchInternalField());

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("k");

    Field<scalar>::operator=
    (
        (alpha_*Tinf_.value() + lambda*TP/delta)
       /(lambda/delta + alpha_ + SMALL)
    );

    fvPatchField<scalar>::evaluate();
}

Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveFluxFvPatchScalarField::
snGrad() const
{
    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);
    scalarField TP(this->patchInternalField());

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("k");

//     scalarField gradient =
//          lambda*alpha_*(Tinf_.value() - TP)
//         /(lambda + alpha_*delta + SMALL);

//     scalarField gradient =
//          alpha_*(Tinf_.value() - TP)
//         /(lambda + alpha_*delta + SMALL);

    return tmp<Field<scalar> >
    (
         alpha_*(Tinf_.value() - TP)
        /(lambda + alpha_*delta + SMALL)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveFluxFvPatchScalarField::
gradientInternalCoeffs() const
{
    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("k");

//     return tmp<Field<scalar> >
//     (
//         -pTraits<scalar>::one*lambda*alpha_/(lambda + alpha_*delta + SMALL)
//     );

    return tmp<Field<scalar> >
    (
        -pTraits<scalar>::one*alpha_/(lambda + alpha_*delta + SMALL)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveFluxFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    scalarField delta(1.0/this->patch().deltaCoeffs() + SMALL);

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("k");

    return alpha_*Tinf_.value()/(lambda + alpha_*delta + SMALL);
}


void Foam::convectiveFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "Tinf", Tinf_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        convectiveFluxFvPatchScalarField
    );
}

// ************************************************************************* //
