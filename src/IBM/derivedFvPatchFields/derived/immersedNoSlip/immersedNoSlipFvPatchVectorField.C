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

#include "immersedNoSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedNoSlipFvPatchVectorField::immersedNoSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    immersedFixedValueFvPatchField<vector>(p, iF)
{
    this->immersedField_ = Zero;
    this->internalValue_ = Zero;
}


Foam::immersedNoSlipFvPatchVectorField::immersedNoSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    immersedFixedValueFvPatchField<vector>(p, iF, dict, false)
{
    this->immersedField_ = Zero;
    this->internalValue_ = Zero;
}


Foam::immersedNoSlipFvPatchVectorField::immersedNoSlipFvPatchVectorField
(
    const immersedNoSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedFixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    this->immersedField_ = Zero;
    this->internalValue_ = Zero;
}


Foam::immersedNoSlipFvPatchVectorField::immersedNoSlipFvPatchVectorField
(
    const immersedNoSlipFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    immersedFixedValueFvPatchField<vector>(ptf, iF)
{
    this->immersedField_ = Zero;
    this->internalValue_ = Zero;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedNoSlipFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    immersedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::immersedNoSlipFvPatchVectorField::setInternal()
{
    if (!this->setInternal_)
    {
        return;
    }
    this->object_.setInternal
    (
        this->internalFieldRef(),
        vector::zero
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        immersedNoSlipFvPatchVectorField
    );
}


// ************************************************************************* //
