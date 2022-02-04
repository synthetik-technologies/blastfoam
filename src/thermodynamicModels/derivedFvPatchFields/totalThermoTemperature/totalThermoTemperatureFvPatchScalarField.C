/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "totalThermoTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalThermoTemperatureFvPatchScalarField::
totalThermoTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    thermoBasePatchField(this->patch(), iF),
    UName_(IOobject::groupName("U", this->phaseName_)),
    phiName_(IOobject::groupName("phi", this->phaseName_)),
    T0_(p.size(), 0.0)
{}


Foam::totalThermoTemperatureFvPatchScalarField::
totalThermoTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    thermoBasePatchField(this->patch(), iF),
    UName_
    (
        dict.lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", this->phaseName_)
        )
    ),
    phiName_
    (
        dict.lookupOrDefault<word>
        (
            "phi",
            IOobject::groupName("phi", this->phaseName_)
        )
    ),
    T0_("T0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(T0_);
    }
}


Foam::totalThermoTemperatureFvPatchScalarField::
totalThermoTemperatureFvPatchScalarField
(
    const totalThermoTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    thermoBasePatchField(this->patch(), iF),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    T0_(mapper(ptf.T0_))
{}


Foam::totalThermoTemperatureFvPatchScalarField::
totalThermoTemperatureFvPatchScalarField
(
    const totalThermoTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    thermoBasePatchField(this->patch(), iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalThermoTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(T0_, T0_);
}


void Foam::totalThermoTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const totalThermoTemperatureFvPatchScalarField& tiptf =
        refCast<const totalThermoTemperatureFvPatchScalarField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
}


void Foam::totalThermoTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalarField psip(this->psi());
    const scalarField gammap(this->gamma());
    const scalarField gM1ByG((gammap - 1.0)/gammap);

    operator==
    (
        T0_/(1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up))
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::totalThermoTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    thermoBasePatchField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "T0", T0_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalThermoTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
