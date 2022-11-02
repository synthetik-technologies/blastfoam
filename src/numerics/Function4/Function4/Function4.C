/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
01-11-2022 Synthetik Applied Technologies : Added Function4
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "Function4.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function4<Type>::Function4(const word& name)
:
    name_(name)
{}


template<class Type>
Foam::Function4<Type>::Function4(const Function4<Type>& de)
:
    tmp<Function4<Type>>::refCount(),
    name_(de.name_)
{}


template<class Type, class Function4Type>
Foam::FieldFunction4<Type, Function4Type>::FieldFunction4
(
    const word& name
)
:
    Function4<Type>(name)
{}


template<class Type, class Function4Type>
Foam::tmp<Foam::Function4<Type>>
Foam::FieldFunction4<Type, Function4Type>::clone() const
{
    return tmp<Function4<Type>>
    (
        new Function4Type(refCast<const Function4Type>(*this))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function4<Type>::~Function4()
{}


template<class Type, class Function4Type>
Foam::FieldFunction4<Type, Function4Type>::~FieldFunction4()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function4<Type>::name() const
{
    return name_;
}


template<class Type, class Function4Type>
Foam::tmp<Foam::Field<Type>> Foam::FieldFunction4<Type, Function4Type>::value
(
    const scalar t,
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
            refCast<const Function4Type>(*this).value
            (
                t,
                x[i],
                y[i],
                z[i]
            );
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function4<Type>::operator=(const Function4<Type>& f)
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
void  Foam::writeEntry(Ostream& os, const Function4<Type>& f4)
{
    writeKeyword(os, f4.name())
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    writeEntry(os, "type", f4.type());

    f4.write(os);

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function4<Type>& f4
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Function4<Type>&)"
    );

    f4.write(os);

    return os;
}


// ************************************************************************* //
