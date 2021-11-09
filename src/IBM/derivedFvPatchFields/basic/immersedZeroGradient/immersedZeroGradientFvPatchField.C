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

#include "immersedZeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedZeroGradientFvPatchField<Type>::immersedZeroGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::immersedZeroGradientFvPatchField<Type>::immersedZeroGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    immersedValueFvPatchField<Type>(p, iF, dict, false)
{
    this->readInternal(dict);
}


template<class Type>
Foam::immersedZeroGradientFvPatchField<Type>::immersedZeroGradientFvPatchField
(
    const immersedZeroGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::immersedZeroGradientFvPatchField<Type>::immersedZeroGradientFvPatchField
(
    const immersedZeroGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(ptf, iF)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedZeroGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->immersedField_ =
        this->object_.patchExternalField(this->internalField());

    immersedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::immersedZeroGradientFvPatchField<Type>::setInternal()
{
    if (this->setInternal_)
    {
        if (this->nSmooth_ > 0)
        {
            this->smoothInternalField(this->internalFieldRef());
        }
        else
        {
            this->setInternalField
            (
                this->internalValue_
            );
        }
    }
}

// ************************************************************************* //
