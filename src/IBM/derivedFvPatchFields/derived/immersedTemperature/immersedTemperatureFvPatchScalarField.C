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

#include "immersedTemperatureFvPatchScalarField.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedTemperatureFvPatchScalarField::immersedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedValueFvPatchField<scalar>(p, iF),
    isEnergy_(iF.member() != "T")
{}


Foam::immersedTemperatureFvPatchScalarField::immersedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    immersedValueFvPatchField<scalar>(p, iF, dict, valueRequired),
    isEnergy_(iF.member() != "T")
{}


Foam::immersedTemperatureFvPatchScalarField::immersedTemperatureFvPatchScalarField
(
    const immersedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    isEnergy_(iF.member() != "T")
{}


Foam::immersedTemperatureFvPatchScalarField::immersedTemperatureFvPatchScalarField
(
    const immersedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedValueFvPatchField<scalar>(ptf, iF),
    isEnergy_(iF.member() != "T")
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::immersedTemperatureFvPatchScalarField::setPatchHE()
{
    if (isEnergy_)
    {
        const labelList cells(this->object_.patchCells());
        const labelListList map(this->object_.patchMap());
        const basicThermo& thermo = basicThermo::lookupThermo(*this);

        scalarField pT(UIndirectList<scalar>(thermo.T(), cells)());
        scalarField he(thermo.he(pT, cells)());

        forAll(map, i)
        {
            forAll(map[i], j)
            {
                this->immersedField_[map[i][j]] = he[i];
            }
        }
    }
}


void Foam::immersedTemperatureFvPatchScalarField::setInternalHE()
{
    if (isEnergy_)
    {
        const labelList cells(this->object_.internalCells());
        const basicThermo& thermo = basicThermo::lookupThermo(*this);
        scalarField T(UIndirectList<scalar>(thermo.T(), cells)());

        this->setInternalField(thermo.he(T, cells));
    }
}
// ************************************************************************* //
