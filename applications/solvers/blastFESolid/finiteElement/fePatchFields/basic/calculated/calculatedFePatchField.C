/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "calculatedFePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::fePatchField<Type>::calculatedType()
{
    return calculatedFePatchField<Type>::typeName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedFePatchField<Type>::calculatedFePatchField
(
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF
)
:
    fePatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedFePatchField<Type>::calculatedFePatchField
(
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF,
    const dictionary& dict
)
:
    fePatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::calculatedFePatchField<Type>::calculatedFePatchField
(
    const calculatedFePatchField<Type>& ptf,
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF,
    const fePatchFieldMapper& mapper
)
:
    fePatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::calculatedFePatchField<Type>::calculatedFePatchField
(
    const calculatedFePatchField<Type>& ptf,
    const DimensionedField<Type, feMesh>& iF
)
:
    fePatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
Foam::autoPtr<Foam::fePatchField<Type>>
Foam::fePatchField<Type>::NewCalculatedType
(
    const fePatchField<Type2>& pf
)
{
    typename fePatchConstructorTable::iterator patchTypeCstrIter =
        fePatchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != fePatchConstructorTablePtr_->end())
    {
        return autoPtr<fePatchField<Type>>
        (
            patchTypeCstrIter()
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
    else
    {
        return autoPtr<fePatchField<Type>>
        (
            new calculatedFePatchField<Type>
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
}


// ************************************************************************* //
