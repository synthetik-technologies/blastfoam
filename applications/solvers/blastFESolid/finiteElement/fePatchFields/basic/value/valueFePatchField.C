/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "valueFePatchField.H"
#include "fePatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::valueFePatchField<Type>::checkFieldSize() const
{
    if (this->size() != this->patch().size())
    {
        FatalErrorInFunction
            << "field does not correspond to patch. " << endl
            << "Field size: " << size() << " patch size: "
            << this->patch().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::valueFePatchField<Type>::valueFePatchField
(
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF
)
:
    fePatchField<Type>(p, iF),
    Field<Type>(p.size())
{}


template<class Type>
Foam::valueFePatchField<Type>::valueFePatchField
(
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fePatchField<Type>(p, iF, dict),
    Field<Type>(p.size())
{
    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else if (!valueRequired)
    {
        Field<Type>::operator=(Zero);
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::valueFePatchField<Type>::valueFePatchField
(
    const valueFePatchField<Type>& ptf,
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF,
    const fePatchFieldMapper& mapper
)
:
    fePatchField<Type>(ptf, p, iF, mapper),
    Field<Type>(mapper(ptf))
{}


template<class Type>
Foam::valueFePatchField<Type>::valueFePatchField
(
    const valueFePatchField<Type>& ptf,
    const DimensionedField<Type, feMesh>& iF
)
:
    fePatchField<Type>(ptf, iF),
    Field<Type>(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::valueFePatchField<Type>::autoMap
(
    const fePatchFieldMapper& m
)
{
    m(*this, *this);
}


template<class Type>
void Foam::valueFePatchField<Type>::rmap
(
    const fePatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap
    (
        refCast<const valueFePatchField<Type>>
        (
            ptf
        ),
        addr
    );
}


template<class Type>
void Foam::valueFePatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInInternalField(iF, *this);

    fePatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::valueFePatchField<Type>::evaluate(const Pstream::commsTypes)
{
    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInInternalField(iF, *this);

    fePatchField<Type>::evaluate();
}


template<class Type>
void Foam::valueFePatchField<Type>::write(Ostream& os) const
{
    fePatchField<Type>::write(os);
    writeEntry(os, "value", static_cast<const Field<Type>&>(*this));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::valueFePatchField<Type>::operator=
(
    const valueFePatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::valueFePatchField<Type>::operator=
(
    const fePatchField<Type>& ptf
)
{
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
void Foam::valueFePatchField<Type>::operator=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::valueFePatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::valueFePatchField<Type>::operator==
(
    const valueFePatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::valueFePatchField<Type>::operator==
(
    const fePatchField<Type>& ptf
)
{
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
void Foam::valueFePatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::valueFePatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// ************************************************************************* //
