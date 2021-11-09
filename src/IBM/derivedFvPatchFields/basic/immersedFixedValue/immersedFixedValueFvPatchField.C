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

#include "immersedFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedFixedValueFvPatchField<Type>::immersedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::immersedFixedValueFvPatchField<Type>::immersedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    immersedValueFvPatchField<Type>(p, iF, dict, valueRequired)
{}


template<class Type>
Foam::immersedFixedValueFvPatchField<Type>::immersedFixedValueFvPatchField
(
    const immersedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::immersedFixedValueFvPatchField<Type>::immersedFixedValueFvPatchField
(
    const immersedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(ptf, iF)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    immersedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::immersedFixedValueFvPatchField<Type>::setInternal()
{
    if (this->setInternal_)
    {
        this->setInternalField(this->internalValue_);
    }
}


template<class Type>
void Foam::immersedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    immersedValueFvPatchField<Type>::write(os);
    writeEntry(os, "immersedValue", this->immersedField_);
}


// ************************************************************************* //
