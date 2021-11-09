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

#include "immersedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedValueFvPatchField<Type>::immersedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedFvPatchField<Type>(p, iF),
    immersedField_(this->object_.nFaces(), Zero)
{}


template<class Type>
Foam::immersedValueFvPatchField<Type>::immersedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    immersedFvPatchField<Type>(p, iF, dict),
    immersedField_(this->object_.nFaces(), Zero)
{
    if (valueRequired)
    {
        if (dict.found("immersedValue"))
        {
            immersedField_ =
                Field<Type>
                (
                    "immersedValue",
                    dict,
                    this->object_.nFaces()
                );
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Essential entry 'immersedValue' missing"
                << exit(FatalIOError);
        }

        if (this->setInternal_ && this->nSmooth_ < 1)
        {
            this->internalValue_ =
                dict.lookup<Type>("internalValue");
        }
    }
}


template<class Type>
Foam::immersedValueFvPatchField<Type>::immersedValueFvPatchField
(
    const immersedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedFvPatchField<Type>(ptf, p, iF, mapper),
    immersedField_(ptf.immersedField_)
{}


template<class Type>
Foam::immersedValueFvPatchField<Type>::immersedValueFvPatchField
(
    const immersedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedFvPatchField<Type>(ptf, iF),
    immersedField_(ptf.immersedField_)
{}


template<class Type>
void Foam::immersedValueFvPatchField<Type>::readInternal
(
    const dictionary& dict
)
{
    if (this->setInternal_ && this->nSmooth_ < 1)
    {
        this->internalValue_ =
            dict.lookup<Type>("internalValue");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper&
)
{}


template<class Type>
void Foam::immersedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& f,
    const labelList&
)
{

    immersedField_ =
        dynamicCast<const immersedValueFvPatchField<Type>&>
        (
            f
        ).immersedField();
}


template<class Type>
void Foam::immersedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (this->setPatchInternal_)
    {
        this->setPatchInternalField(immersedField_);
    }
    immersedFvPatchField<Type>::updateCoeffs();
    this->setInternal();
}


template<class Type>
void
Foam::immersedValueFvPatchField<Type>::addForcing
(
    Field<Type>& F,
    const Field<scalar>& alphaRho,
    const Field<Type>& old,
    const Field<Type>& RHS,
    const scalar& dt
) const
{
    tmp<Field<Type>> interpF
    (
        (
            immersedField_*this->object_.interpolateTo(alphaRho)
          - this->object_.interpolateTo(old)
        )/dt
      + this->object_.interpolateTo(RHS)

    );
    this->object_.interpolateFrom(interpF(), F);
}


template<class Type>
void Foam::immersedValueFvPatchField<Type>::write(Ostream& os) const
{
    immersedFvPatchField<Type>::write(os);
}

// ************************************************************************* //
