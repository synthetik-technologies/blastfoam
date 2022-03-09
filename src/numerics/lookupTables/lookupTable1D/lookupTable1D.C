/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "lookupTable1D.H"
#include "tableReader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D()
:
    modType_("none"),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    index_(0),
    f_(0.0)
{
    setMod(modType_, modFunc_, invModFunc_);
    setMod(modXType_, modXFunc_, invModXFunc_);
    setInterp(interpType_, interpFunc_);
}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D
(
    const dictionary& dict,
    const word& xName,
    const word& name,
    bool canRead
)
:
    modType_("none"),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    index_(0),
    f_(0.0)
{
    read(dict, xName, name, canRead);
}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D
(
    const Field<scalar>& x,
    const Field<Type>& data,
    const word& xMod,
    const word& mod,
    const word& interpolationScheme,
    const bool isReal
)
:
    modType_(mod),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_(xMod),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_(interpolationScheme),
    interpFunc_(nullptr),
    xValues_(x),
    xModValues_(x),
    data_(data),
    index_(0),
    f_(0.0)
{
    set(x, data, modXType_, modType_, interpType_, isReal);
}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D
(
    const Field<scalar>& x,
    const word& xMod,
    const word& interpolationScheme,
    const bool isReal
)
:
    modType_("none"),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_(xMod),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_(interpolationScheme),
    interpFunc_(nullptr),
    xValues_(),
    xModValues_(),
    data_(),
    index_(0),
    f_(0.0)
{
    setX(x, modXType_, interpType_, isReal);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable1D<Type>::~lookupTable1D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable1D<Type>::set
(
    const Field<scalar>& x,
    const Field<Type>& data,
    const bool isReal
)
{
    setX(x, isReal);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::set
(
    const Field<scalar>& x,
    const Field<Type>& data,
    const word& xMod,
    const word& mod,
    const word& interpolationScheme,
    const bool isReal
)
{
    setX(x, xMod, interpolationScheme, isReal);
    setData(data, mod, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::setX
(
    const Field<scalar>& x,
    const bool isReal
)
{
    if (isReal)
    {
        xValues_ = x;
        xModValues_.resize(x.size());
        forAll(x, i)
        {
            xModValues_[i] = modXFunc_(x[i]);
        }
    }
    else
    {
        xModValues_ = x;
        xValues_.resize(x.size());
        forAll(x, i)
        {
            xValues_[i] = invModXFunc_(x[i]);
        }
    }
}


template<class Type>
void Foam::lookupTable1D<Type>::setX
(
    const Field<scalar>& x,
    const word& xMod,
    const word& interpolationScheme,
    const bool isReal
)
{
    modXType_ = xMod;
    setMod(modXType_, modXFunc_, invModXFunc_);
    setX(x, isReal);

    interpType_ = interpolationScheme;
    setInterp(interpType_, interpFunc_);
}


template<class Type>
void Foam::lookupTable1D<Type>::setData
(
    const Field<Type>& data,
    const bool isReal
)
{
    data_ = data;
    if (isReal)
    {
        forAll(data_, i)
        {
            data_[i] = modFunc_(data[i]);
        }
    }
}


template<class Type>
void Foam::lookupTable1D<Type>::setData
(
    const Field<Type>& data,
    const word& mod,
    const bool isReal
)
{
    modType_ = mod;
    setMod(modType_, modFunc_, invModFunc_);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::update(const scalar x) const
{
    if (x <= xValues_[0])
    {
        index_ = 0;
        f_ = 0.0;
        return;
    }
    else if (x >= xValues_.last())
    {
        index_ = xValues_.size() - 2;
        f_ = 1.0;
        return;
    }
    for (index_ = 1; index_ < xModValues_.size(); index_++)
    {
        if (x <= xValues_[index_])
        {
            index_--;
            break;
        }
    }
    f_ =
        linearWeight
        (
            modXFunc_(x),
            xModValues_[index_],
            xModValues_[index_ + 1]
        );
    return;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::lookupTable1D<Type>::realData() const
{
#ifdef FULL_DEBUG
    if (!invModFunc_)
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    tmp<Field<Type>> tmpF(new Field<Type>(data_));
    Field<Type>& f = tmpF.ref();
    forAll(f, i)
    {
        f[i] = invModFunc_(f[i]);
    }
    return tmpF;
}


template<class Type>
Type Foam::lookupTable1D<Type>::lookup(const scalar x) const
{
#ifdef FULL_DEBUG
    if (!invModFunc_)
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    update(x);
    return invModFunc_(interpFunc_(modXFunc_(x), index_, xModValues_, data_));
}


template<class Type>
Type Foam::lookupTable1D<Type>::dFdX(const scalar x) const
{
#ifdef FULL_DEBUG
    if (!invModFunc_)
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    update(x);

    scalar fm(data_[index_]);
    scalar fp(data_[index_ + 1]);

    return
        (invModFunc_(fp) - invModFunc_(fm))
       /(xValues_[index_ + 1] - xValues_[index_]);
}


template<class Type>
Type Foam::lookupTable1D<Type>::d2FdX2(const scalar x) const
{
#ifdef FULL_DEBUG
    if (!invModFunc_)
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    update(modXFunc_(x));
    if (index_ == 0)
    {
        index_++;
    }

    scalar ym(invModFunc_(data_[index_-1]));
    scalar yi(invModFunc_(data_[index_]));
    scalar yp(invModFunc_(data_[index_+1]));

    const scalar& xm(xValues_[index_-1]);
    const scalar& xi(xValues_[index_]);
    const scalar& xp(xValues_[index_+1]);

    return
        ((yp - yi)/(xp - xi) - (yi - ym)/(xi - xm))/(xp - xm);
}


template<class Type>
void Foam::lookupTable1D<Type>::read
(
    const dictionary& dict,
    const word& xName,
    const word& name,
    const bool canRead
)
{
    interpType_ =
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp");
    setInterp(interpType_, interpFunc_);

    List<List<string>> table;
    if (dict.found("file"))
    {
        table = read2DTable
        (
            dict.lookup<fileName>("file"),
            dict.lookupOrDefault<string>("delim", ","),
            dict.lookupOrDefault<label>("startRow", 0),
            dict.lookupOrDefault<Switch>("flipTable", false)
        );
    }

    readComponent<scalar>
    (
        dict,
        xName,
        table,
        modXType_,
        xValues_,
        xModValues_,
        modXFunc_,
        invModXFunc_,
        canRead
    );

    Field<Type> data;
    readComponent<Type>
    (
        dict,
        name,
        table,
        modType_,
        data,
        data_,
        modFunc_,
        invModFunc_,
        canRead
    );
}

// ************************************************************************* //
