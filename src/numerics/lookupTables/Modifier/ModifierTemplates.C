/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "Modifier.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Modifier<Type>> Foam::Modifier<Type>::New
(
    const word& modifierType
)
{
    DebugInfo
        << "Selecting " << pTraits<Type>::typeName
        << " modifier: " << modifierType << endl;

    typename nullConstructorTable::iterator cstrIter =
        nullConstructorTablePtr_->find(modifierType);

    if (cstrIter == nullConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << Modifier<Type>::typeName << " type "
            << modifierType << nl << nl
            << "Valid " << Modifier<Type>::typeName << " types are : " << endl
            << nullConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<Modifier<Type>>(cstrIter()());
}


template<class Type>
Foam::autoPtr<Foam::Modifier<Type>> Foam::Modifier<Type>::New
(
    const word& modifierType,
    const dictionary& dict
)
{
    DebugInfo
        << "Selecting " << pTraits<Type>::typeName
        << " modifier: " << modifierType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modifierType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << Modifier<Type>::typeName << " type "
            << modifierType << nl << nl
            << "Valid " << Modifier<Type>::typeName << " types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<Modifier<Type>>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Modifier<Type>::Mod(UList<Type>& f) const
{
    forAll(f, i)
    {
        f[i] = this->operator()(f[i]);
    }
}


template<class Type>
void Foam::Modifier<Type>::Inv(UList<Type>& f) const
{
    forAll(f, i)
    {
        f[i] = this->inv(f[i]);
    }
}

// ************************************************************************* //
