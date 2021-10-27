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

#include "lookupTable2D.H"
#include "tableReader.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

template<class Type>
Type Foam::lookupTable2D<Type>::getValue
(
    const label ij,
    const scalar f,
    const List<Type>& xy
) const
{
    if (ij >= xy.size())
    {
        return xy.last();
    }

    return xy[ij] + f*(xy[ij+1] - xy[ij]);
}


template<class Type>
bool Foam::lookupTable2D<Type>::checkUniform(const List<scalar>& xy) const
{
    scalar dxy = xy[1] - xy[0];

    for (label i = 2; i < xy.size(); i++)
    {
        if (mag(xy[i] - xy[i-1] - dxy) > small)
        {
            return false;
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D()
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    xValues_(),
    yValues_(),
    i_(0),
    j_(0),
    fx_(0),
    fy_(0)
{}


template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    xValues_(),
    yValues_(),
    i_(0),
    j_(0),
    fx_(0),
    fy_(0)
{
    read(dict, xName, yName, name);
}


template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<Field<Type>>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const word& interpolationScheme,
    const bool isReal
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    interpFunc_(nullptr),
    data_(data),
    xModValues_(x),
    yModValues_(y),
    xValues_(x),
    yValues_(y),
    i_(0),
    j_(0),
    fx_(0),
    fy_(0)
{
    set(x, y, data, modXType, modYType, modType, interpolationScheme, isReal);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable2D<Type>::~lookupTable2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable2D<Type>::set
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<Field<Type>>& data,
    const bool isReal
)
{
    setX(x, isReal);
    setY(y, isReal);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::set
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<Field<Type>>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const word& interpolationScheme,
    const bool isReal
)
{
    setX(x, modXType, isReal);
    setY(y, modYType, isReal);
    setData(data, modType, isReal);
    setInterp(interpolationScheme, interpFunc_);
}


template<class Type>
void Foam::lookupTable2D<Type>::setX
(
    const Field<scalar>& x,
    const word& modXType,
    const bool isReal
)
{
    setMod(modXType, modXFunc_, invModXFunc_);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::setX
(
    const Field<scalar>& x,
    const bool isReal
)
{
    if (!isReal)
    {
        forAll(xValues_, j)
        {
            xValues_[j] = invModXFunc_(xModValues_[j]);
        }
    }
    else
    {
        forAll(xValues_, j)
        {
            xModValues_[j] = modXFunc_(xValues_[j]);
        }
    }

    if (checkUniform(xModValues_))
    {
        findXIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        findXIndex_ = &lookupTable2D::findNonuniformIndexes;
    }
}


template<class Type>
void Foam::lookupTable2D<Type>::setY
(
    const Field<scalar>& x,
    const word& modXType,
    const bool isReal
)
{
    setMod(modXType, modXFunc_, invModXFunc_);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::setY
(
    const Field<scalar>& y,
    const bool isReal
)
{
    if (!isReal)
    {
        forAll(yValues_, j)
        {
            yValues_[j] = invModYFunc_(yModValues_[j]);
        }
    }
    else
    {
        forAll(yValues_, j)
        {
            yModValues_[j] = modYFunc_(yValues_[j]);
        }
    }

    if (checkUniform(yModValues_))
    {
        findYIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        findYIndex_ = &lookupTable2D::findNonuniformIndexes;
    }
}


template<class Type>
void Foam::lookupTable2D<Type>::setData
(
    const Field<Field<Type>>& data,
    const bool isReal
)
{
    data_ = data;

    if (isReal)
    {
        forAll(data_, i)
        {
            forAll(data_[i], j)
            {
                data_[i][j] = modFunc_(data_[i][j]);
            }
        }
    }
}


template<class Type>
void Foam::lookupTable2D<Type>::setData
(
    const Field<Field<Type>>& data,
    const word& modType,
    const bool isReal
)
{
    setMod(modType, modFunc_, invModFunc_);
    setData(data, isReal);
}


template<class Type>
Foam::tmp<Foam::Field<Foam::Field<Type>>>
Foam::lookupTable2D<Type>::realData() const
{
    tmp<Field<Field<Type>>> tmpf(new Field<Field<Type>>(data_));
    Field<Field<Type>>& f = tmpf.ref();
    forAll(f, i)
    {
        forAll(f[i], j)
        {
            f[i][j] = invModFunc_(f[i][j]);
        }
    }
    return tmpf;
}


template<class Type>
void Foam::lookupTable2D<Type>::update(const scalar x, const scalar y) const
{
    updateX(x);
    updateY(y);
}


template<class Type>
void Foam::lookupTable2D<Type>::updateX(const scalar x) const
{
    i_ = findXIndex_(modXFunc_(x), xModValues_);
}


template<class Type>
void Foam::lookupTable2D<Type>::updateY(const scalar y) const
{
    j_ = findYIndex_(modYFunc_(y), yModValues_);

}


template<class Type>
void Foam::lookupTable2D<Type>::updateXWeight(const scalar x) const
{
    fx_ = linearWeight(modXFunc_(x), xModValues_[i_], xModValues_[i_ + 1]);
}


template<class Type>
void Foam::lookupTable2D<Type>::updateYWeight(const scalar y) const
{
    fy_ = linearWeight(modYFunc_(y), yModValues_[j_], yModValues_[j_ + 1]);

}


template<class Type>
Type Foam::lookupTable2D<Type>::lookup(const scalar x, const scalar y) const
{
    update(x, y);
    return
        invModFunc_
        (
            interpFunc_
            (
                modXFunc_(x), modYFunc_(y),
                i_, j_,
                xModValues_, yModValues_,
                data_
            )
        );
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdX(const scalar x, const scalar y) const
{
    update(x, y);
    updateYWeight(y);

    return
        (
            invModFunc_
            (
                data_[i_+1][j_]*(1.0 - fy_)
              + data_[i_+1][j_+1]*fy_
            )
          - invModFunc_
            (
                data_[i_][j_]*(1.0 - fy_)
              + data_[i_][j_+1]*fy_
            )
        )/(xValues_[i_+1] - xValues_[i_]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdY(const scalar x, const scalar y) const
{
    update(x, y);
    updateXWeight(x);
    return
        (
            invModFunc_
            (
                data_[i_][j_+1]*(1.0 - fx_)
              + data_[i_+1][j_+1]*fx_
            )
          - invModFunc_
            (
                data_[i_][j_]*(1.0 - fx_)
              + data_[i_+1][j_]*fx_
            )
        )/(yValues_[j_+1] - yValues_[j_]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdX2(const scalar x, const scalar y) const
{
    update(x, y);
    updateYWeight(y);

    if (i_ == 0)
    {
        i_++;
    }

    const Type gmm(invModFunc_(data_[i_-1][j_]));
    const Type gm(invModFunc_(data_[i_][j_]));
    const Type gpm(invModFunc_(data_[i_+1][j_]));

    const Type gmp(invModFunc_(data_[i_-1][j_+1]));
    const Type gp(invModFunc_(data_[i_][j_+1]));
    const Type gpp(invModFunc_(data_[i_+1][j_+1]));

    const scalar xm(xValues_[i_-1]);
    const scalar xi(xValues_[i_]);
    const scalar xp(xValues_[i_+1]);

    const Type gPrimepm((gpm - gm)/(xp - xi));
    const Type gPrimemm((gm - gmm)/(xi - xm));
    const Type gPrimepp((gpp - gp)/(xp - xi));
    const Type gPrimemp((gp - gmp)/(xi - xm));

    return
        (1.0 - fy_)*(gPrimepm - gPrimemm)/(0.5*(xp - xm))
      + fy_*(gPrimepp - gPrimemp)/(0.5*(xp - xm));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdY2(const scalar x, const scalar y) const
{
    update(x, y);
    updateXWeight(x);

    if (j_ == 0)
    {
        j_++;
    }

    const Type gmm(invModFunc_(data_[i_][j_-1]));
    const Type gm(invModFunc_(data_[i_][j_]));
    const Type gmp(invModFunc_(data_[i_][j_+1]));

    const Type gpm(invModFunc_(data_[i_+1][j_-1]));
    const Type gp(invModFunc_(data_[i_+1][j_]));
    const Type gpp(invModFunc_(data_[i_+1][j_+1]));

    const scalar ym(yValues_[j_-1]);
    const scalar yi(yValues_[j_]);
    const scalar yp(yValues_[j_+1]);

    const Type gPrimemp((gmp - gm)/(yp - yi));
    const Type gPrimemm((gm - gmm)/(yi - ym));
    const Type gPrimepp((gpp - gp)/(yp - yi));
    const Type gPrimepm((gp - gpm)/(yi - ym));

    return
        (1.0 - fx_)*(gPrimemp - gPrimemm)/(0.5*(yp - ym))
      + fx_*(gPrimepp - gPrimepm)/(0.5*(yp - ym));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdXdY(const scalar x, const scalar y) const
{
    update(x, y);

    const Type gmm(invModFunc_(data_[i_][j_]));
    const Type gmp(invModFunc_(data_[i_][j_+1]));
    const Type gpm(invModFunc_(data_[i_+1][j_]));
    const Type gpp(invModFunc_(data_[i_+1][j_+1]));

    const scalar xm(xValues_[i_]);
    const scalar xp(xValues_[i_+1]);

    const scalar ym(yValues_[j_]);
    const scalar yp(yValues_[j_+1]);

    return ((gpp - gmp)/(xp - xm) - (gpm - gmm)/(xp - xm))/(yp - ym);
}


template<class Type>
void Foam::lookupTable2D<Type>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name
)
{
    word interpolationScheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );
    setInterp(interpolationScheme, interpFunc_);

    readComponent
    (
        dict,
        xName,
        xValues_,
        xModValues_,
        modXFunc_,
        invModXFunc_,
        findXIndex_
    );
    readComponent
    (
        dict,
        yName,
        yValues_,
        yModValues_,
        modYFunc_,
        invModYFunc_,
        findYIndex_
    );

    data_.resize(xValues_.size());
    forAll(data_, i)
    {
        data_[i].resize(yValues_.size());
    }

    const dictionary& fDict(dict.subDict(name + "Coeffs"));
    setMod
    (
        fDict.lookupOrDefault<word>("mod", "none"),
        modFunc_,
        invModFunc_
    );
    Switch isReal(fDict.lookupOrDefault<Switch>("isReal", true));

    fileName file(fDict.lookup<word>("file"));
    read2DTable
    (
        file,
        fDict.lookupOrDefault<string>("delim", ","),
        data_,
        fDict.lookupOrDefault<Switch>("flipTable", true)
    );
    if (isReal)
    {
        forAll(xValues_, i)
        {
            forAll(yValues_, j)
            {
                data_[i][j] = modFunc_(data_[i][j]);
            }
        }
    }
}


template<class Type>
void Foam::lookupTable2D<Type>::readComponent
(
    const dictionary& parentDict,
    const word& name,
    Field<scalar>& values,
    Field<scalar>& modValues,
    modFuncType& modFunc,
    modFuncType& invModFunc,
    findIndexFunc& findIndex
)
{
    Switch isReal = true;
    if (parentDict.found(name))
    {
        values = parentDict.lookup<Field<Type>>(name);
        modValues.resize(values.size());
        isReal =
            parentDict.lookupOrDefault<Switch>(name + "isReal", true);
        setMod
        (
            parentDict.lookupOrDefault<word>(name + "Mod", "none"),
            modFunc,
            invModFunc
        );
    }
    else if (parentDict.found(name + "Coeffs"))
    {
        const dictionary& dict(parentDict.subDict(name + "Coeffs"));
        isReal = dict.lookupOrDefault<Switch>("isReal", true);
        setMod
        (
            dict.lookupOrDefault<word>("mod", "none"),
            modFunc,
            invModFunc
        );

        if (dict.found("file"))
        {
            fileName file(dict.lookup("file"));
            read1DTable
            (
                file,
                dict.lookupOrDefault<string>("delim", ","),
                values
            );
        }
        else
        {
            label ny = dict.lookup<label>("n");
            Type dy = dict.lookup<Type>("delta");
            Type miny = dict.lookup<Type>("min");

            values.resize(ny);
            forAll(values, j)
            {
                values[j] = miny + dy*j;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Either a list of values or " << name << "Coeffs must" << nl
            << "be provided for " << name << endl
            << abort(FatalError);
    }

    if (!isReal)
    {
        modValues = values;
        forAll(values, i)
        {
            values[i] = invModFunc(values[i]);
        }
    }
    else
    {
        modValues.resize(values.size());
        forAll(values, i)
        {
            modValues[i] = modFunc(values[i]);
        }
    }
    if (checkUniform(modValues))
    {
        findIndex = &findUniformIndexes;
    }
    else
    {
        findIndex = &findNonuniformIndexes;
    }
}


// ************************************************************************* //
