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

#include "List2D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::List2D<Type>::List2D(const label m, const label n)
:
    List<Type>(m*n),
    m_(m),
    n_(n),
    own_(true)
{}


template<class Type>
Foam::List2D<Type>::List2D(const label m, const label n, const Type& t)
:
    List<Type>(m*n, t),
    m_(m),
    n_(n),
    own_(true)
{}


template<class Type>
template<class Type2>
Foam::List2D<Type>::List2D(const label m, const label n, const Type2& t)
:
    List<Type>(m*n, t),
    m_(m),
    n_(n),
    own_(true)
{}


template<class Type>
Foam::List2D<Type>::List2D(Type* v, const label m, const label n)
:
    List<Type>(),
    m_(m),
    n_(n),
    own_(false)
{
    UList<Type> lst(v, m_*n_);
    UList<Type>::shallowCopy(lst);
}


template<class Type>
Foam::List2D<Type>::List2D(const List2D<Type>& lst)
:
    List<Type>(lst),
    m_(lst.m_),
    n_(lst.n_),
    own_(true)
{}


template<class Type>
template<class ListListType>
Foam::List2D<Type>::List2D(const ListListType& lst)
:
    List<Type>(),
    m_(),
    n_(),
    own_(true)
{
    operator=(lst);
}


template<class Type>
Foam::List2D<Type>::List2D(List2D<Type>&& lst)
:
    List<Type>(),
    m_(),
    n_(),
    own_(lst.own_)
{
    if (!own_)
    {
        m_ = lst.m_;
        n_ = lst.n_;
        UList<Type>::shallowCopy(lst);
    }
    else
    {
        transfer(lst);
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::List2D<Type>::~List2D()
{
    if (!own_)
    {
        UList<Type> lst(nullptr, 0);
        UList<Type>::shallowCopy(lst);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::List2D<Type>::resize(const label m, const label n)
{
    m_ = m;
    n_ = n;
    const label newSize = m*n;
    if (List<Type>::size() < newSize)
    {
        List<Type>::setSize(newSize);
    }
}


template<class Type>
void Foam::List2D<Type>::resize(const label m, const label n, const Type& t)
{
    m_ = m;
    n_ = n;
    const label newSize = m*n;
    if (List<Type>::size() < newSize)
    {
        List<Type>::setSize(newSize, t);
    }
}


template<class Type>
void Foam::List2D<Type>::setSize(const label m, const label n)
{
    resize(m, n);
}


template<class Type>
void Foam::List2D<Type>::setSize(const label m, const label n, const Type& t)
{
    resize(m, n, t);
}


template<class Type>
void Foam::List2D<Type>::transfer(List2D<Type>& lst)
{
    m_ = lst.m_;
    n_= lst.n_;
    List<Type>::transfer(lst);
}


template<class Type>
void Foam::List2D<Type>::flip()
{
    List2D<Type> newLst(n_, m_);
    for (label i = 0; i < m_; i++)
    {
        for (label j = 0; j < n_; j++)
        {
            newLst(j, i) = operator()(i, j);
        }
    }
    operator=(move(newLst));
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::List2D<Type>::operator=(const List2D<Type>& lst)
{
    m_ = lst.m_;
    n_ = lst.n_;
    static_cast<List<Type>&>(*this) = lst;
}


template<class Type>
template<class Type2>
void Foam::List2D<Type>::operator=(const List2D<Type2>& lst)
{
    m_ = lst.m_;
    n_ = lst.n_;
    static_cast<List<Type>&>(*this) = lst;
}


template<class Type>
template<class ListListType>
void Foam::List2D<Type>::operator=(const ListListType& ll)
{
    resize(ll.size(), ll.size() ? ll[0].size() : 0);
    forAll(ll, i)
    {
        if (n_ != ll[i].size())
        {
            FatalErrorInFunction
                << "inconsistent row sizes" << endl
                << exit(FatalError);
        }
        forAll(ll[i], j)
        {
            operator()(i, j) = ll[i][j];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List2DIO.C"

// ************************************************************************* //
