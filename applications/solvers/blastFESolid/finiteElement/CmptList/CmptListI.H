/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "CmptList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::CmptList<T>::CmptList(UList<T>& lst, const label cmpt)
:
    UList<T>(lst.begin(), lst.size()),
    cmpt_(cmpt)
{}


template<class T>
Foam::CmptList<T>::CmptList(const UList<T>& lst, const label cmpt)
:
    UList<T>(const_cast<T*>(lst.cbegin()), lst.size()),
    cmpt_(cmpt)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::scalar& Foam::CmptList<T>::first() const
{
    return component(UList<T>::first(), this->cmpt_);
}


template<class T>
inline Foam::scalar& Foam::CmptList<T>::first()
{
    return setComponent(UList<T>::first(), this->cmpt_);
}


template<class T>
inline const Foam::scalar& Foam::CmptList<T>::last() const
{
    return component(UList<T>::last(), this->cmpt_);
}


template<class T>
inline Foam::scalar& Foam::CmptList<T>::last()
{
    return setComponent(UList<T>::last(), cmpt_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const Foam::scalar& Foam::CmptList<T>::operator[](const label i) const
{
    return component(UList<T>::operator[](i), this->cmpt_);
}


template<class T>
inline Foam::scalar& Foam::CmptList<T>::operator[](const label i)
{
    return setComponent(UList<T>::operator[](i), this->cmpt_);
}


template<class T>
inline void Foam::CmptList<T>::operator=(const scalar& v)
{
    forAll((*this), i)
    {
        operator[](i) = v;
    }
}


template<class T>
inline void Foam::CmptList<T>::operator=(const zero)
{
    forAll((*this), i)
    {
        operator[](i) = Zero;
    }
}


// ************************************************************************* //
