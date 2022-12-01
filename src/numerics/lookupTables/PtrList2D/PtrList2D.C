/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "PtrList2D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PtrList2D<Type>::PtrList2D(const label m, const label n)
:
    PtrList<Type>(m*n),
    m_(m),
    n_(n),
    own_(true)
{}


template<class Type>
Foam::PtrList2D<Type>::PtrList2D(Type** v, const label m, const label n)
:
    PtrList<Type>(m*n),
    m_(m),
    n_(n),
    own_(false)
{
    PtrList<Type>& lst(*this);
    forAll(lst, i)
    {
        lst(i) = v[i];
    }
}


template<class Type>
Foam::PtrList2D<Type>::PtrList2D(PtrList2D<Type>&& lst)
:
    PtrList<Type>(lst.m_*lst.n_),
    m_(),
    n_(),
    own_(lst.own_)
{
    if (!own_)
    {
        m_ = lst.m_;
        n_ = lst.n_;
        PtrList<Type>& ptrlst(lst);
        forAll(lst, i)
        {
            PtrList<Type>::set(i, &ptrlst[i]);
        }
    }
    else
    {
        transfer(lst);
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::PtrList2D<Type>::~PtrList2D()
{
    if (!own_)
    {
        UPtrList<Type>::clear();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PtrList2D<Type>::resize(const label m, const label n)
{
    m_ = m;
    n_ = n;
    const label newSize = m*n;
    if (PtrList<Type>::size() < newSize)
    {
        PtrList<Type>::setSize(newSize);
    }
}


template<class Type>
void Foam::PtrList2D<Type>::setSize(const label m, const label n)
{
    resize(m, n);
}


template<class Type>
void Foam::PtrList2D<Type>::transfer(PtrList2D<Type>& lst)
{
    m_ = lst.m_;
    n_= lst.n_;
    PtrList<Type>::transfer(lst);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::PtrList2D<Type>::operator=(const PtrList2D<Type>& lst)
{
    m_ = lst.m_;
    n_ = lst.n_;
    static_cast<PtrList<Type>&>(*this) = lst;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrList2DIO.C"

// ************************************************************************* //
