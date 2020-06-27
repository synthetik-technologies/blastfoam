/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Changed class name
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "blastMixedEnergyCalculatedTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        blastMixedEnergyCalculatedTemperatureFvPatchScalarField,
        0
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::
blastMixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    heRefValue_(p.size()),
    heRefGrad_(p.size()),
    heValueFraction_(p.size())
{}


Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::
blastMixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF),
    heRefValue_(p.size()),
    heRefGrad_(p.size()),
    heValueFraction_(p.size())
{
    calculatedFvPatchScalarField::evaluate();
}


Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::
blastMixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const blastMixedEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchScalarField(ptf, p, iF, mapper),
    heRefValue_(mapper(ptf.heRefValue_)),
    heRefGrad_(mapper(ptf.heRefGrad_)),
    heValueFraction_(mapper(ptf.heValueFraction_))
{}


Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::
blastMixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const blastMixedEnergyCalculatedTemperatureFvPatchScalarField& ptf
)
:
    calculatedFvPatchScalarField(ptf),
    heRefValue_(ptf.heRefValue_),
    heRefGrad_(ptf.heRefGrad_),
    heValueFraction_(ptf.heValueFraction_)
{}


Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::
blastMixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const blastMixedEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(ptf, iF),
    heRefValue_(ptf.heRefValue_),
    heRefGrad_(ptf.heRefGrad_),
    heValueFraction_(ptf.heValueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    calculatedFvPatchScalarField::autoMap(m);
    m(heRefValue_, heRefValue_);
    m(heRefGrad_, heRefGrad_);
    m(heValueFraction_, heValueFraction_);
}


void Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    calculatedFvPatchScalarField::rmap(ptf, addr);

    const blastMixedEnergyCalculatedTemperatureFvPatchScalarField& mptf =
        refCast<const blastMixedEnergyCalculatedTemperatureFvPatchScalarField>(ptf);

    heRefValue_.rmap(mptf.heRefValue_, addr);
    heRefGrad_.rmap(mptf.heRefGrad_, addr);
    heValueFraction_.rmap(mptf.heValueFraction_, addr);
}


// ************************************************************************* //
