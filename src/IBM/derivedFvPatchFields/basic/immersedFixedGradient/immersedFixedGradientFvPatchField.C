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

#include "immersedFixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedFixedGradientFvPatchField<Type>::immersedFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(p, iF),
    gradient_(this->object_.nFaces(), Zero)
{}


template<class Type>
Foam::immersedFixedGradientFvPatchField<Type>::immersedFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool gradientRequired
)
:
    immersedValueFvPatchField<Type>(p, iF, dict, false),
    gradient_(this->object_.nFaces(), Zero)
{
    this->readInternal(dict);
    if (gradientRequired)
    {
        if (dict.found("gradient"))
        {
            gradient_ =
                Field<Type>
                (
                    "gradient",
                    dict,
                    gradient_.size()
                );
            evaluate();
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Essential entry 'gradient' missing"
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::immersedFixedGradientFvPatchField<Type>::immersedFixedGradientFvPatchField
(
    const immersedFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedValueFvPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_)
{}


template<class Type>
Foam::immersedFixedGradientFvPatchField<Type>::immersedFixedGradientFvPatchField
(
    const immersedFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedFixedGradientFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper&
)
{}


template<class Type>
void Foam::immersedFixedGradientFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& f,
    const labelList& l
)
{
    immersedValueFvPatchField<Type>::rmap(f, l);
    gradient_ =
        dynamicCast<const immersedFixedGradientFvPatchField<Type>&>
        (
            f
        ).gradient();
}


template<class Type>
void Foam::immersedFixedGradientFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    this->immersedField_ =
        this->object_.patchExternalField(this->internalField())
      + gradient_/this->object_.deltaCoeffs();

    immersedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::immersedFixedGradientFvPatchField<Type>::setInternal()
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


template<class Type>
void Foam::immersedFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    immersedValueFvPatchField<Type>::write(os);
    writeEntry(os, "gradient", gradient_);
}


// ************************************************************************* //
