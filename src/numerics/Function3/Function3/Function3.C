/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "Function3.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function3<Type>::Function3(const word& name)
:
    name_(name)
{}


template<class Type>
Foam::Function3<Type>::Function3(const Function3<Type>& de)
:
    tmp<Function3<Type>>::refCount(),
    name_(de.name_)
{}


template<class Type, class Function3Type>
Foam::FieldFunction3<Type, Function3Type>::FieldFunction3
(
    const word& name
)
:
    Function3<Type>(name)
{}


template<class Type, class Function3Type>
Foam::tmp<Foam::Function3<Type>>
Foam::FieldFunction3<Type, Function3Type>::clone() const
{
    return tmp<Function3<Type>>
    (
        new Function3Type(refCast<const Function3Type>(*this))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function3<Type>::~Function3()
{}


template<class Type, class Function3Type>
Foam::FieldFunction3<Type, Function3Type>::~FieldFunction3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function3<Type>::name() const
{
    return name_;
}


template<class Type, class Function3Type>
Foam::tmp<Foam::Field<Type>> Foam::FieldFunction3<Type, Function3Type>::value
(
    const scalarField& x,
    const scalarField& y,
    const scalarField& z
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] =
            refCast<const Function3Type>(*this).value
            (
                x[i],
                y[i],
                z[i]
            );
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function3<Type>::operator=(const Function3<Type>& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class Type>
void  Foam::writeEntry(Ostream& os, const Function3<Type>& f3)
{
    writeKeyword(os, f3.name())
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    writeEntry(os, "type", f3.type());

    f3.write(os);

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function3<Type>& f3
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Function3<Type>&)"
    );

    f3.write(os);

    return os;
}


// ************************************************************************* //
