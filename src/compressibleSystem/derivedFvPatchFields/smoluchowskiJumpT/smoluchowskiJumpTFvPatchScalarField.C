/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "smoluchowskiJumpTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluidThermoModel.H"
#include "mathematicalConstants.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Pr_(1.0),
    accommodationCoeff_(1.0),
    Twall_(p.size(), 0.0),
    gamma_(1.4)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const smoluchowskiJumpTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Pr_(ptf.Pr_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    Twall_(ptf.Twall_),
    gamma_(ptf.gamma_)
{}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Pr_(dict.lookupType<scalar>("Pr")),
    accommodationCoeff_(readScalar(dict.lookup("accommodationCoeff"))),
    Twall_("Twall", dict, p.size()),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.4))
{
    if
    (
        mag(accommodationCoeff_) < small
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "unphysical accommodationCoeff specified"
            << "(0 < accommodationCoeff <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const smoluchowskiJumpTFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    Pr_(ptpsf.Pr_),
    accommodationCoeff_(ptpsf.accommodationCoeff_),
    Twall_(ptpsf.Twall_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::smoluchowskiJumpTFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::smoluchowskiJumpTFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void Foam::smoluchowskiJumpTFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fluidThermoModel& thermo =
        db().lookupObject<fluidThermoModel>
        (
            IOobject::groupName("basicThermo", internalField().group())
        );
    const label patchi = patch().index();
    const scalarField& pmu = thermo.mu(patchi);
    const scalarField& prho = thermo.rho().boundaryField()[patchi];
    const scalarField& pT = thermo.T().boundaryField()[patchi];

    Field<scalar> C2
    (
        pmu/prho
       *sqrt
        (
            Foam::constant::thermodynamic::RR/thermo.W(patchi)*pT
           *constant::mathematical::piByTwo
        )
       *2.0*gamma_/Pr_/(gamma_ + 1.0)
       *(2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C2));
    refValue() = Twall_;
    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::smoluchowskiJumpTFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "accommodationCoeff", accommodationCoeff_);
    writeEntry(os, "Pr", Pr_);
    writeEntry(os, "Twall", Twall_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        smoluchowskiJumpTFvPatchScalarField
    );
}


// ************************************************************************* //
