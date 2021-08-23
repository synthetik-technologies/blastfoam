/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
#include "DynamicList.H"
#include "Field.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

void Foam::lookupTable1D::readTable
(
    const fileName& file,
    const label xi,
    const label yi,
    const bool isReal
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

    scalarField xTmp;
    scalarField xModTmp;
    scalarField fTmp;
    string line;
    while (is.good())
    {
        is.getLine(line);
        if (line[0] == '#')
        {
            continue;
        }
        line.replaceAll(delim_, " ");
        line = '(' + line + ')';

        IStringStream isLine(line);

        scalarList lineVals(isLine);
        if (lineVals.size())
        {
            if (isReal)
            {
                xTmp.append(lineVals[xi]);
                xModTmp.append(modXFunc_(lineVals[xi]));
                fTmp.append(modFunc_(lineVals[yi]));
            }
            else
            {
                xModTmp.append(lineVals[xi]);
                xTmp.append(invModXFunc_(lineVals[xi]));
                fTmp.append(lineVals[yi]);
            }
        }
    }
    xValues_ = xTmp;
    xModValues_ = xModTmp;
    data_ = fTmp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupTable1D::lookupTable1D()
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpFunc_(nullptr),
    index_(0),
    f_(0.0)
{}


Foam::lookupTable1D::lookupTable1D
(
    const dictionary& dict,
    const word& xName,
    const word& yName
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
    word mod(dict.lookupOrDefault<word>("mod", "none"));
    word xMod(dict.lookupOrDefault<word>(xName + "Mod", "none"));
    word interpolationScheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );
    Switch isReal(dict.lookupOrDefault<Switch>("isReal", true));

    setMod(mod, modFunc_, invModFunc_);
    setMod(xMod, modXFunc_, invModXFunc_);
    setInterp(interpolationScheme, interpFunc_);

    if (dict.found("delim"))
    {
        delim_ = dict.lookup<string>("delim");
    }

    if (dict.found("file"))
    {
        fileName file(dict.lookup("file"));
        label xi(dict.lookupOrDefault<label>(xName + "Col", 0));
        label yi(dict.lookupOrDefault<label>(yName + "Col", 1));
        readTable(file, xi, yi, isReal);
    }
    else
    {
        xValues_ = dict.lookup<scalarField>(xName);
        xModValues_ = xValues_;
        data_ = dict.lookup<scalarField>(yName);
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
}


Foam::lookupTable1D::lookupTable1D
(
    const fileName& file,
    const word& mod,
    const word& xMod,
    const word& interpolationScheme
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
    setMod(mod, modFunc_, invModFunc_);
    setMod(xMod, modXFunc_, invModXFunc_);
    setInterp(interpolationScheme, interpFunc_);
    readTable(file);
}


Foam::lookupTable1D::lookupTable1D
(
    const scalarField& x,
    const scalarField& data,
    const word& mod,
    const word& xMod,
    const word& interpolationScheme,
    const bool correct
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
    setMod(mod, modFunc_, invModFunc_);
    setMod(xMod, modXFunc_, invModXFunc_);
    setInterp(interpolationScheme, interpFunc_);
    if (!correct)
    {
        forAll(x, i)
        {
            xModValues_[i] = invModXFunc_(x[i]);
        }
    }
    else
    {
        forAll(x, i)
        {
            xValues_[i] = modXFunc_(x[i]);
            data_[i] = modFunc_(x[i]);
        }
    }
}


Foam::lookupTable1D::lookupTable1D
(
    const scalarField& x,
    const word& xMod,
    const word& interpolationScheme,
    const bool correct
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    interpFunc_(nullptr),
    xValues_(x),
    xModValues_(x),
    data_(),
    index_(0),
    f_(0.0)
{
    setMod(xMod, modXFunc_, invModXFunc_);
    setInterp(interpolationScheme, interpFunc_);
    if (!correct)
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
            xValues_[i] = invModXFunc_(x[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lookupTable1D::~lookupTable1D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lookupTable1D::update(const scalar& x) const
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
    f_ = interpFunc_(modXFunc_(x), xModValues_[index_], xModValues_[index_ + 1]);
    return;
}


Foam::scalar Foam::lookupTable1D::lookup(const scalar& x) const
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
    return invModFunc_(data_[index_] + f_*(data_[index_+1] - data_[index_]));
}


Foam::scalar
Foam::lookupTable1D::reverseLookup(const scalar& yin) const
{
#ifdef FULL_DEBUG
    if (!modFunc_)
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    scalar y(modFunc_(yin));
    if (y < data_[0])
    {
        index_ = 0;
    }
    if (y > data_.last())
    {
        index_ = data_.size() - 2;
    }
    else
    {
        for (index_ = 0; index_ < data_.size(); index_++)
        {
            if (y < data_[index_])
            {
                index_--;
                break;
            }
        }
    }

    const scalar& ym(data_[index_]);
    const scalar& yp(data_[index_+1]);

    f_ = interpFunc_(y, ym, yp);

    return invModXFunc_
    (
        xModValues_[index_]
      + f_*(xModValues_[index_+1] - xModValues_[index_])
    );
}


Foam::scalar Foam::lookupTable1D::dFdX(const scalar& x) const
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


Foam::scalar Foam::lookupTable1D::d2FdX2(const scalar& x) const
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


void Foam::lookupTable1D::set
(
    const scalarField& x,
    const scalarField& data,
    const word& mod,
    const word& xMod,
    const word& interpolationScheme,
    const bool inReal
)
{
    setMod(mod, modFunc_, invModFunc_);
    setMod(xMod, modXFunc_, invModXFunc_);
    setInterp(interpolationScheme, interpFunc_);
    if (inReal)
    {
        xValues_ = x;
        xModValues_.resize(x.size());
        data_ = data;
        forAll(x, i)
        {
            xModValues_[i] = modXFunc_(x[i]);
        }
    }
    else
    {
        xModValues_ = x;
        xValues_.resize(x.size());
        data_ = data;
        forAll(x, i)
        {
            xValues_[i] = invModXFunc_(x[i]);
            data_[i] = invModFunc_(data[i]);
        }
    }
}

void Foam::lookupTable1D::setX
(
    const scalarField& x,
    const bool inReal
)
{
    data_.clear();

    if (inReal)
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

// ************************************************************************* //
