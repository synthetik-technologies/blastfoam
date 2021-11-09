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

#include "immersedMixedFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedMixedFvPatchField<Type>::immersedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(p, iF),
    refValue_(this->object_.nFaces(), Zero),
    refGrad_(this->object_.nFaces(), Zero),
    valueFraction_(this->object_.nFaces(), Zero)
{}


template<class Type>
Foam::immersedMixedFvPatchField<Type>::immersedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    immersedValueFvPatchField<Type>(p, iF, dict),
    refValue_("refValue", dict, this->object_.nFaces()),
    refGrad_("refGrad", dict, this->object_.nFaces()),
    valueFraction_("valueFraction", dict, this->object_.nFaces())
{
    evaluate();
}


template<class Type>
Foam::immersedMixedFvPatchField<Type>::immersedMixedFvPatchField
(
    const immersedMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedValueFvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


template<class Type>
Foam::immersedMixedFvPatchField<Type>::immersedMixedFvPatchField
(
    const immersedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedValueFvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedMixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper&
)
{}


template<class Type>
void Foam::immersedMixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& f,
    const labelList& l
)
{
    immersedValueFvPatchField<Type>::rmap(f, l);
    const immersedMixedFvPatchField<Type>& pf =
        dynamicCast<const immersedMixedFvPatchField<Type>&>(f);
    refValue_ = pf.refValue();
    refGrad_ = pf.refGrad();
    valueFraction_ = pf.valueFraction();
}


template<class Type>
void Foam::immersedMixedFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    this->immersedField_ =
        valueFraction_*refValue_
      + (1.0 - valueFraction_)
       *(
            this->object_.patchExternalField(this->internalField())
          + refGrad_/this->object_.deltaCoeffs()
        );

    immersedValueFvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::immersedMixedFvPatchField<Type>::setInternal()
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
void Foam::immersedMixedFvPatchField<Type>::write(Ostream& os) const
{
    immersedValueFvPatchField<Type>::write(os);
    writeEntry(os, "refValue", refValue_);
    writeEntry(os, "refGrad", refGrad_);
    writeEntry(os, "valueFraction", valueFraction_);
}

// ************************************************************************* //
