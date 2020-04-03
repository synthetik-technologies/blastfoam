/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "blastGradientEnergyCalculatedTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        blastGradientEnergyCalculatedTemperatureFvPatchScalarField,
        0
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::
blastGradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    heGradient_(p.size())
{}


Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::
blastGradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF),
    heGradient_(p.size())
{
    calculatedFvPatchScalarField::evaluate();
}


Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::
blastGradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const blastGradientEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchScalarField(ptf, p, iF, mapper),
    heGradient_(mapper(ptf.heGradient_))
{}


Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::
blastGradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const blastGradientEnergyCalculatedTemperatureFvPatchScalarField& ptf
)
:
    calculatedFvPatchScalarField(ptf),
    heGradient_(ptf.heGradient_)
{}


Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::
blastGradientEnergyCalculatedTemperatureFvPatchScalarField
(
    const blastGradientEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(ptf, iF),
    heGradient_(ptf.heGradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    calculatedFvPatchScalarField::autoMap(m);
    m(heGradient_, heGradient_);
}


void Foam::blastGradientEnergyCalculatedTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    calculatedFvPatchScalarField::rmap(ptf, addr);

    const blastGradientEnergyCalculatedTemperatureFvPatchScalarField& mptf =
        refCast<const blastGradientEnergyCalculatedTemperatureFvPatchScalarField>
        (
            ptf
        );

    heGradient_.rmap(mptf.heGradient_, addr);
}


// ************************************************************************* //
