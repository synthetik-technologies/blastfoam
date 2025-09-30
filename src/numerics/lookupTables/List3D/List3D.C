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

#include "List3D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::List3D<Type>::List3D(const label m, const label n, const label l)
:
    List<Type>(m*n*l),
    m_(m),
    n_(n),
    l_(l),
    own_(true)
{}


template<class Type>
Foam::List3D<Type>::List3D
(
    const label m,
    const label n,
    const label l,
    const Type& t
)
:
    List<Type>(m*n*l, t),
    m_(m),
    n_(n),
    l_(l),
    own_(true)
{}


template<class Type>
template<class Type2>
Foam::List3D<Type>::List3D
(
    const label m,
    const label n,
    const label l,
    const Type2& t
)
:
    List<Type>(m*n*l, t),
    m_(m),
    n_(n),
    l_(l),
    own_(true)
{}


template<class Type>
Foam::List3D<Type>::List3D(Type* v, const label m, const label n, const label l)
:
    List<Type>(),
    m_(m),
    n_(n),
    l_(l),
    own_(false)
{
    UList<Type> lst(v, m_*n_*l_);
    UList<Type>::shallowCopy(lst);
}


template<class Type>
Foam::List3D<Type>::List3D(const List3D& lst)
:
    List<Type>(lst),
    m_(lst.m_),
    n_(lst.n_),
    l_(lst.l_),
    own_(true)
{}


template<class Type>
Foam::List3D<Type>::List3D(List3D<Type>&& lst)
:
    List<Type>(),
    m_(),
    n_(),
    l_(),
    own_(lst.own_)
{
    if (!own_)
    {
        UList<Type>::shallowCopy(lst);
    }
    else
    {
        transfer(lst);
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::List3D<Type>::~List3D()
{
    if (!own_)
    {
        UList<Type> lst(nullptr, 0);
        UList<Type>::shallowCopy(lst);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::List3D<Type>::resize(const label m, const label n, const label l)
{
    m_ = m;
    n_ = n;
    l_ = l;
    const label newSize = m*n*l;
    if (List<Type>::size() < newSize)
    {
        List<Type>::setSize(newSize);
    }
}


template<class Type>
void Foam::List3D<Type>::resize
(
    const label m,
    const label n,
    const label l,
    const Type& t
)
{
    m_ = m;
    n_ = n;
    l_ = l;
    const label newSize = m*n*l;
    if (List<Type>::size() < newSize)
    {
        List<Type>::setSize(newSize, t);
    }
}


template<class Type>
void Foam::List3D<Type>::setSize(const label m, const label n, const label l)
{
    resize(m, n, l);
}


template<class Type>
void Foam::List3D<Type>::setSize
(
    const label m,
    const label n,
    const label l,
    const Type& t
)
{
    resize(m, n, l, t);
}


template<class Type>
void Foam::List3D<Type>::transfer(List3D<Type>& lst)
{
    m_ = lst.m_;
    n_= lst.n_;
    l_ = lst.l_;
    own_ = lst.own_;
    if (!lst.own_)
    {
        UList<Type>::shallowCopy(lst);
    }
    else
    {
        List<Type>::transfer(lst);
    }
}





template<class Type>
void Foam::List3D<Type>::flip()
{
    List3D<Type> newLst(l_, n_, m_);
    for (label k = 0; k < l_; k++)
    {
        for (label j = 0; j < n_; j++)
        {
            for (label i = 0; i < m_; i++)
            {
                newLst[k][j][i] = operator()(i, j, k);
            }
        }
    }
    operator=(newLst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::List3D<Type>::operator=(const List3D<Type>& lst)
{
    m_ = lst.m_;
    n_ = lst.n_;
    l_ = lst.l_;
    own_ = lst.own_;
    static_cast<List<Type>&>(*this) = lst;
}


template<class Type>
template<class Type2>
void Foam::List3D<Type>::operator=(const List3D<Type2>& lst)
{
    m_ = lst.m_;
    n_ = lst.n_;
    l_ = lst.l_;
    static_cast<List<Type>&>(*this) = lst;
}


template<class Type>
template<class ListListListType>
void Foam::List3D<Type>::operator=(const ListListListType& ll)
{
    label s1 = ll.size();
    label s2 = s1 ? ll[0].size() : 0;
    label s3 = s2 ? ll[0][0].size() : 0;
    resize(s1, s2, s3);
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
            forAll(ll[i][j], k)
            {
                operator()(i, j, k) = ll[i][j][k];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List3DIO.C"

// ************************************************************************* //
