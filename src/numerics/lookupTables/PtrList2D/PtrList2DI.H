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
inline bool Foam::PtrList2D<Type>::checkIIndex(const label i) const
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
inline bool Foam::PtrList2D<Type>::checkJIndex(const label j) const
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
inline Foam::PtrList2D<Type>::PtrList2D()
:
    m_(0),
    n_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::label Foam::PtrList2D<Type>::index
(
    const label i,
    const label j
) const
{
    return i*n_ + j;
}

template<class Type>
inline Foam::label Foam::PtrList2D<Type>::m() const
{
    return m_;
}


template<class Type>
inline Foam::label Foam::PtrList2D<Type>::n() const
{
    return n_;
}


template<class Type>
inline Foam::label Foam::PtrList2D<Type>::size() const
{
    return m_*n_;
}


template<class Type>
inline Foam::autoPtr<Type> Foam::PtrList2D<Type>::set
(
    const labelVector2D& ij,
    Type* v
)
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return PtrList<Type>::set(index(ij.x(), ij.y()), v);
}


template<class Type>
inline Foam::autoPtr<Type> Foam::PtrList2D<Type>::set
(
    const label i,
    const label j,
    Type* v
)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return PtrList<Type>::set(index(i, j), v);
}


template<class Type>
inline Foam::autoPtr<Type> Foam::PtrList2D<Type>::set
(
    const labelVector2D& ij,
    autoPtr<Type>& v
)
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return PtrList<Type>::set(index(ij.x(), ij.y()), v);
}


template<class Type>
inline Foam::autoPtr<Type> Foam::PtrList2D<Type>::set
(
    const label i,
    const label j,
    autoPtr<Type>& v
)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return PtrList<Type>::set(index(i, j), v);
}


template<class Type>
inline Foam::autoPtr<Type> Foam::PtrList2D<Type>::set
(
    const labelVector2D& ij,
    tmp<Type>& v
)
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return PtrList<Type>::set(index(ij.x(), ij.y()), v);
}


template<class Type>
inline Foam::autoPtr<Type> Foam::PtrList2D<Type>::set
(
    const label i,
    const label j,
    tmp<Type>& v
)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return PtrList<Type>::set(index(i, j), v);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::PtrList2D<Type>::operator()
(
    const labelVector2D& ij
) const
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return PtrList<Type>::operator[](index(ij.x(), ij.y()));
}


template<class Type>
inline Type& Foam::PtrList2D<Type>::operator()
(
    const labelVector2D& ij
)
{
    #ifdef FULLDEBUG
    checkIIndex(ij.x());
    checkJIndex(ij.y());
    #endif
    return PtrList<Type>::operator[](index(ij.x(), ij.y()));
}


template<class Type>
inline const Type& Foam::PtrList2D<Type>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return PtrList<Type>::operator[](index(i, j));
}


template<class Type>
inline Type& Foam::PtrList2D<Type>::operator()(const label i, const label j)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    checkJIndex(j);
    #endif
    return PtrList<Type>::operator[](index(i, j));
}


template<class Type>
inline Foam::UPtrList<Type> Foam::PtrList2D<Type>::operator()(const label i) const
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    #endif
    UPtrList<Type> row(n_);
    const PtrList<Type>& lst(*this);
    PtrList<Type>& ptrlst(const_cast<PtrList<Type>&>(lst));
    forAll(row, rowi)
    {
        row.set(rowi, &ptrlst[i*n_ + rowi]);
    }
    return row;
}


template<class Type>
inline Foam::UPtrList<Type> Foam::PtrList2D<Type>::operator()(const label i)
{
    #ifdef FULLDEBUG
    checkIIndex(i);
    #endif
    UPtrList<Type> row(n_);
    forAll(row, rowi)
    {
        row.set(rowi, &PtrList<Type>::operator[](i*n_ + rowi));
    }
    return row;
}


template<class Type>
inline Foam::UPtrList<Type> Foam::PtrList2D<Type>::operator[](const label i) const
{
    return operator()(i);
}


template<class Type>
inline Foam::UPtrList<Type> Foam::PtrList2D<Type>::operator[](const label i)
{
    return operator()(i);
}


// ************************************************************************* //
