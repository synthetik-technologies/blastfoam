/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "dynamicThermoPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicThermoPressureFvPatchScalarField, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicThermoPressureFvPatchScalarField::dynamicThermoPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    thermoBasePatchField(p, iF),
    p0_(p.size(), 0)
{}


Foam::dynamicThermoPressureFvPatchScalarField::dynamicThermoPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    thermoBasePatchField(p, iF, dict),
    p0_("p0", dict, p.size())
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
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::dynamicThermoPressureFvPatchScalarField::dynamicThermoPressureFvPatchScalarField
(
    const dynamicThermoPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    thermoBasePatchField(p, iF),
    p0_(mapper(ptf.p0_))
{}


Foam::dynamicThermoPressureFvPatchScalarField::dynamicThermoPressureFvPatchScalarField
(
    const dynamicThermoPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    thermoBasePatchField(tppsf.patch(), iF),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicThermoPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(p0_, p0_);
}


void Foam::dynamicThermoPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const dynamicThermoPressureFvPatchScalarField& tiptf =
        refCast<const dynamicThermoPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::dynamicThermoPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const scalarField& Kp
)
{
    if (updated())
    {
        return;
    }


    // High-speed compressible flow
    const scalarField psip(this->psi());
    const scalarField gammap(this->gamma());
    const scalarField gM1ByG((gammap - 1.0)/gammap);
    operator==(p0p/pow(scalar(1) + psip*gM1ByG*Kp, 1/gM1ByG));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::dynamicThermoPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    thermoBasePatchField::write(os);
    writeEntry(os, "p0", p0_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
