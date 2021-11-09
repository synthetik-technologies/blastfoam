/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "immersedFixedTemperatureFvPatchScalarField.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedFixedTemperatureFvPatchScalarField::immersedFixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedTemperatureFvPatchScalarField(p, iF)
{}


Foam::immersedFixedTemperatureFvPatchScalarField::immersedFixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedTemperatureFvPatchScalarField(p, iF, dict)
{
    if (this->setInternal_)
    {
        this->internalValue_ = dict.lookup<scalar>("internalValue");
    }
}


Foam::immersedFixedTemperatureFvPatchScalarField::immersedFixedTemperatureFvPatchScalarField
(
    const immersedFixedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedTemperatureFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::immersedFixedTemperatureFvPatchScalarField::immersedFixedTemperatureFvPatchScalarField
(
    const immersedFixedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedTemperatureFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedFixedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (!this->updated())
    {
        return;
    }
    if (this->isEnergy_)
    {
        const basicThermo& thermo = basicThermo::lookupThermo(*this);
        const label patchi = this->patch().index();
        const immersedValueFvPatchField& pT
        (
            dynamicCast<const immersedValueFvPatchField>
            (
                thermo.T().boundaryField()[patchi]
            )
        );
        const_cast<immersedValueFvPatchField&>(pT).evaluate();
        this->setPatchHE();
    }
}


void Foam::immersedFixedTemperatureFvPatchScalarField::setInternal()
{
    if (!this->setInternal_)
    {
        return;
    }
    if (this->isEnergy_)
    {
        this->setInternalHE();
    }
    else
    {
        this->setInternalField(this->internalValue_);
    }
}


void Foam::immersedFixedTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    immersedTemperatureFvPatchScalarField::write(os);
    writeEntry(os, "immersedValue", this->immersedField_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        immersedFixedTemperatureFvPatchScalarField
    );
}


// ************************************************************************* //
