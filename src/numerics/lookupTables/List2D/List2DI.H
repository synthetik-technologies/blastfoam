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


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class Type>
inline bool Foam::List2D<Type>::checkIIndex(const label i) const
{
    if (i < 0 || i > m_)
    {
        FatalErrorInFunction
            << "First index " << i << " out of range 0 ... " << m_
            << abort(FatalError);
        return false;
    }
    return true;
}


template<class Type>
inline bool Foam::List2D<Type>::checkJIndex(const label j) const
{
    if (j < 0 || j > n_)
    {
        FatalErrorInFunction
            << "Second index " << j << " out of range 0 ... " << n_
            << abort(FatalError);
        return false;
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::List2D<Type>::List2D()
:
    m_(0),
    n_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::label Foam::List2D<Type>::index
(
    const label i,
    const label j
) const
{
    return i*n_ + j;
}

template<class Type>
inline Foam::label Foam::List2D<Type>::m() const
{
    return m_;
}


template<class Type>
inline Foam::label Foam::List2D<Type>::n() const
{
    return n_;
}


template<class Type>
inline Foam::label Foam::List2D<Type>::size() const
{
    return m_*n_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::List2D<Type>::operator()
(
    const labelVector2D& ij
) const
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return List<Type>::operator[](index(ij.x(), ij.y()));
}


template<class Type>
inline Type& Foam::List2D<Type>::operator()
(
    const labelVector2D& ij
)
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return List<Type>::operator[](index(ij.x(), ij.y()));
}


template<class Type>
inline const Type& Foam::List2D<Type>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return List<Type>::operator[](index(i, j));
}


template<class Type>
inline Type& Foam::List2D<Type>::operator()(const label i, const label j)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return List<Type>::operator[](index(i, j));
}


template<class Type>
inline Foam::UList<Type> Foam::List2D<Type>::operator()(const label i) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    #endif
    return UList<Type>(const_cast<Type*>(this->cdata())+i*n_, n_);
}


template<class Type>
inline Foam::UList<Type> Foam::List2D<Type>::operator()(const label i)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    #endif
    return UList<Type>(this->data()+i*n_, n_);
}


template<class Type>
inline Foam::UList<Type> Foam::List2D<Type>::operator[](const label i) const
{
    return operator()(i);
}


template<class Type>
inline Foam::UList<Type> Foam::List2D<Type>::operator[](const label i)
{
    return operator()(i);
}


// ************************************************************************* //
