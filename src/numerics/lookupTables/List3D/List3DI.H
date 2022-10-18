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
inline bool Foam::List3D<Type>::checkIIndex(const label i) const
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
inline bool Foam::List3D<Type>::checkJIndex(const label j) const
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


template<class Type>
inline bool Foam::List3D<Type>::checkKIndex(const label k) const
{
    if (k < 0 || k > l_)
    {
        FatalErrorInFunction
            << "Third index " << k << " out of range 0 ... " << l_
            << abort(FatalError);
        return false;
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::List3D<Type>::List3D()
:
    m_(0),
    n_(0),
    l_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::label Foam::List3D<Type>::index
(
    const label i,
    const label j,
    const label k
) const
{
    return (i*n_ + j)*l_ + k;
}

template<class Type>
inline Foam::label Foam::List3D<Type>::m() const
{
    return m_;
}


template<class Type>
inline Foam::label Foam::List3D<Type>::n() const
{
    return n_;
}


template<class Type>
inline Foam::label Foam::List3D<Type>::l() const
{
    return l_;
}


template<class Type>
inline Foam::label Foam::List3D<Type>::size() const
{
    return m_*n_*l_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::List3D<Type>::operator()(const labelVector& ijk) const
{
    #ifdef FULLDEBUG
    checkIIndex(ijk.x());
    checkJIndex(ijk.y());
    checkKIndex(ijk.z());
    #endif
    return List<Type>::operator[](index(ijk.x(), ijk.y(), ijk.z()));
}


template<class Type>
inline Type& Foam::List3D<Type>::operator()(const labelVector& ijk)
{
    #ifdef FULLDEBUG
    checkIIndex(ijk.x());
    checkJIndex(ijk.y());
    checkKIndex(ijk.z());
    #endif
    return List<Type>::operator[](index(ijk.x(), ijk.y(), ijk.z()));
}

template<class Type>
inline const Type& Foam::List3D<Type>::operator()
(
    const label i,
    const label j,
    const label k
) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    checkKIndex(k);
    #endif
    return List<Type>::operator[](index(i, j, k));
}


template<class Type>
inline Type& Foam::List3D<Type>::operator()
(
    const label i,
    const label j,
    const label k
)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    checkKIndex(k);
    #endif
    return List<Type>::operator[](index(i, j, k));
}


template<class Type>
inline Foam::UList<Type> Foam::List3D<Type>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return UList<Type>(const_cast<Type*>(this->cdata())+index(i, j, 0), l_);
}


template<class Type>
inline Foam::UList<Type> Foam::List3D<Type>::operator()
(
    const label i,
    const label j
)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return UList<Type>(const_cast<Type*>(this->cdata())+index(i, j, 0), l_);
}


template<class Type>
inline Foam::List2D<Type> Foam::List3D<Type>::operator()(const label i) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    #endif
    return move(List2D<Type>(const_cast<Type*>(this->cdata())+i*n_*l_, n_, l_));
}


template<class Type>
inline Foam::List2D<Type> Foam::List3D<Type>::operator()(const label i)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    #endif
    return move(List2D<Type>(this->data()+i*n_*l_, n_, l_));
}


template<class Type>
inline Foam::List2D<Type> Foam::List3D<Type>::operator[](const label i) const
{
    return operator()(i);
}


template<class Type>
inline Foam::List2D<Type> Foam::List3D<Type>::operator[](const label i)
{
    return operator()(i);
}


// ************************************************************************* //
