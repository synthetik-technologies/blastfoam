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
    isEnergy_(iF.group() != "T")
{}


Foam::immersedTemperatureFvPatchScalarField::immersedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedValueFvPatchField<scalar>(p, iF),
    isEnergy_(iF.group() != "T")
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
    isEnergy_(iF.group() != "T")
{}


Foam::immersedTemperatureFvPatchScalarField::immersedTemperatureFvPatchScalarField
(
    const immersedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedValueFvPatchField<scalar>(ptf, iF),
    isEnergy_(iF.group() != "T")
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (isEnergy_)
    {
        const basicThermo& thermo =
            this->internalField().mesh().lookupObject<basicThermo>
            (
                IOobject::groupName
                (
                    basicThermo::dictName,
                    this->internalField().group()
                )
            );

        const labelList pE(this->object_.patchExternalCells());
        const labelList pI(this->object_.patchInternalCells());
        labelList cells(pE.size());
        Field<scalar> T(this->immersedField_);
        labelList map(pE.size(), -1);

        label I = 0;
        forAll(pE, i)
        {
            label celli = pE[i];
            if (celli >= 0)
            {
                map[I] = i;
                T[I] = T[i];
                cells[I] = celli;
                I++;
            }
            else
            {
                celli = pI[i];
                if (celli >= 0)
                {
                    map[I] = i;
                    T[I] = T[i];
                    cells[I] = celli;
                    I++;
                }
            }
        }
        if (I == 0)
        {
            return;
        }
        map.resize(I);
        cells.resize(I);
        T.resize(I);
        scalarField he(thermo.he(T, cells));
        this->immersedField_ = sum(he)/scalar(I);
        forAll(map, i)
        {
            this->immersedField_[map[i]] = he[i];
        }
    }

    immersedValueFvPatchField<scalar>::updateCoeffs();
    setInternal();

    if (!this->setInternal_)
    {
        return;
    }
}


// ************************************************************************* //
