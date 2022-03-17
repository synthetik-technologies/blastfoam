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
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D()
:
    modType_("none"),
    needMod_(false),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    needXMod_(false),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0.0)
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
    needMod_(false),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    needXMod_(false),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0.0)
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
    needMod_(false),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_(xMod),
    needXMod_(false),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_(interpolationScheme),
    interpFunc_(nullptr),
    data_(data),
    xModValues_(x),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0.0)
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
    needMod_(false),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_(xMod),
    needXMod_(false),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpType_(interpolationScheme),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0.0)
{
    setX(x, modXType_, interpType_, isReal);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable1D<Type>::~lookupTable1D()
{
    if (xValuesPtr_ != &xModValues_)
    {
        deleteDemandDrivenData(xValuesPtr_);
    }
    if (realDataPtr_ != &data_)
    {
        deleteDemandDrivenData(realDataPtr_);
    }
}


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
    setX(x, xMod, isReal);
    setData(data, mod, isReal);
    setInterp(interpolationScheme, interpFunc_);

}


template<class Type>
void Foam::lookupTable1D<Type>::setX
(
    const Field<scalar>& x,
    const bool isReal
)
{
    if (xValuesPtr_ != &xModValues_ && xValuesPtr_ != nullptr)
    {
        deleteDemandDrivenData(xValuesPtr_);
    }

    if (!needXMod_)
    {
        xModValues_ = x;
        xValuesPtr_ = &xModValues_;
        return;
    }

    xValuesPtr_ = new scalarField(x);
    xModValues_ = x;


    if (isReal)
    {
        forAll(x, i)
        {
            xModValues_[i] = modXFunc_(x[i]);
        }
    }
    else
    {
        forAll(x, i)
        {
            (*xValuesPtr_)[i] = invModXFunc_(x[i]);
        }
    }
}


template<class Type>
void Foam::lookupTable1D<Type>::setX
(
    const Field<scalar>& x,
    const word& xMod,
    const bool isReal
)
{
    modXType_ = xMod;
    needXMod_ = xMod != "none";

    setMod(modXType_, modXFunc_, invModXFunc_);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::setData
(
    const Field<Type>& data,
    const bool isReal
)
{
    if (realDataPtr_ != &data_ && realDataPtr_ != nullptr)
    {
        deleteDemandDrivenData(realDataPtr_);
    }

    if (!needMod_)
    {
        data_ = data;
        realDataPtr_ = &data_;
        return;
    }

    realDataPtr_ = new Field<Type>(data);
    data_ = data;
    if (isReal)
    {
        forAll(data, i)
        {
            data_[i] = modFunc_(data[i]);
        }
    }
    else
    {
        forAll(data, i)
        {
            (*realDataPtr_)[i] = modFunc_(data[i]);
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
    needMod_ = mod != "none";

    setMod(modType_, modFunc_, invModFunc_);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::update(const scalar x) const
{
    index_ = findLower(xValues(), x);
    interpFunc_
    (
        modXFunc_(x),
        index_,
        xModValues_,
        indices_,
        weights_
    );
    return;
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

    Type modf = weights_[0]*data_[indices_[0]];
    for (label i = 1; i < indices_.size(); i++)
    {
        modf += weights_[i]*data_[indices_[i]];
    }
    return invModFunc_(modf);
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
       /(xValues()[index_ + 1] - xValues()[index_]);
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

    const scalar& xm(xValues()[index_-1]);
    const scalar& xi(xValues()[index_]);
    const scalar& xp(xValues()[index_+1]);

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

    scalarField x;
    bool isReal = readComponent<scalar>
    (
        dict,
        xName,
        modXType_,
        x,
        table
    );
    setX(x, modXType_, isReal);

    Field<Type> data;
    isReal = readComponent<Type>
    (
        dict,
        name,
        modType_,
        data,
        table
    );
    setData(data, modType_, isReal);
}

// ************************************************************************* //
