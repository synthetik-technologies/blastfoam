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

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable1D<Type>::readTable
(
    const fileName& file,
    const string& delim,
    const label xi,
    const label yi
)
{
    fileName fNameExpanded(file);
    fNameExpanded.expand();

    // Open a stream and check it
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fNameExpanded));
    ISstream& is = isPtr();
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file" << file << nl
            << exit(FatalIOError);
    }

    scalarList xTmp;
    scalarList xModTmp;
    scalarList fTmp;
    string line;
    while (is.good())
    {
        is.getLine(line);
        if (line[0] == '#')
        {
            continue;
        }
        line.replaceAll(delim, " ");
        line = '(' + line + ')';

        IStringStream isLine(line);

        scalarList lineVals(isLine);
        if (lineVals.size())
        {
            xTmp.append(lineVals[xi]);
            fTmp.append(lineVals[yi]);
        }
    }
    xValues_ = xTmp;
    data_ = fTmp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D()
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpFunc_(nullptr),
    index_(0),
    f_(0.0)
{}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D
(
    const dictionary& dict,
    const word& xName,
    const word& name
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpFunc_(nullptr),
    index_(0),
    f_(0.0)
{
    read(dict, xName, name);
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
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpFunc_(nullptr),
    xValues_(x),
    xModValues_(x),
    data_(data),
    index_(0),
    f_(0.0)
{
    set(x, data, xMod, mod, interpolationScheme, isReal);
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
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpFunc_(nullptr),
    xValues_(),
    xModValues_(),
    data_(),
    index_(0),
    f_(0.0)
{
    setX(x, xMod, interpolationScheme, isReal);
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
    setMod(xMod, modXFunc_, invModXFunc_);
    setX(x, isReal);
    setInterp(interpolationScheme, interpFunc_);
}


template<class Type>
void Foam::lookupTable1D<Type>::setData
(
    const Field<Type>& data,
    const bool isReal
)
{
    data_ = data;

    if (!isReal)
    {
        data_ = data;
        forAll(data_, i)
        {
            data_[i] = invModFunc_(data[i]);
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
    setMod(mod, modFunc_, invModFunc_);
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
    f_ = linearWeight(modXFunc_(x), xModValues_[index_], xModValues_[index_ + 1]);
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
        f[i] - invModFunc_(f[i]);
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
    const word& name
)
{
    word mod(dict.lookupOrDefault<word>(name + "Mod", "none"));
    word xMod(dict.lookupOrDefault<word>(xName + "Mod", "none"));
    word interpolationScheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );
    Switch isReal(dict.lookupOrDefault<Switch>("isReal", true));

    setMod(mod, modFunc_, invModFunc_);
    setMod(xMod, modXFunc_, invModXFunc_);
    setInterp(interpolationScheme, interpFunc_);

    if (dict.found("file"))
    {
        readTable
        (
            dict.lookup<fileName>("file"),
            dict.lookupOrDefault<string>("delim", ","),
            dict.lookupOrDefault<label>(xName + "Col", 0),
            dict.lookupOrDefault<label>(name + "Col", 1)
        );
    }
    else
    {
        xValues_ = dict.lookup<scalarList>(xName);
        data_ = dict.lookup<scalarList>(name);

    }

    xModValues_ = xValues_;
    if (!isReal)
    {
        forAll(xValues_, i)
        {
            xValues_[i] = invModXFunc_(xValues_[i]);
        }
    }
    else
    {
        forAll(xValues_, i)
        {
            xModValues_[i] = modXFunc_(xValues_[i]);
            data_[i] = modFunc_(data_[i]);
        }
    }
}

// ************************************************************************* //
