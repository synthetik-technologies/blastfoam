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
#include "demandDrivenData.H"

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
    modType_("none"),
    needMod_(false),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    needXMod_(false),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYType_("none"),
    needYMod_(false),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    ij_(0, 0),
    indices_(0),
    weights_(0)
{}


template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name,
    const bool canRead
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
    modYType_("none"),
    needYMod_(false),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    ij_(0, 0),
    indices_(0),
    weights_(0)
{
    read(dict, xName, yName, name, canRead);
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
    modType_("none"),
    needMod_(false),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    needXMod_(false),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYType_("none"),
    needYMod_(false),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
    data_(data),
    xModValues_(x),
    yModValues_(y),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    ij_(0, 0),
    indices_(0),
    weights_(0)
{
    set(x, y, data, modXType, modYType, modType, interpolationScheme, isReal);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable2D<Type>::~lookupTable2D()
{
    if (xValuesPtr_ != &xModValues_)
    {
        deleteDemandDrivenData(xValuesPtr_);
    }
    if (yValuesPtr_ != &yModValues_)
    {
        deleteDemandDrivenData(yValuesPtr_);
    }
    if (realDataPtr_ != &data_)
    {
        deleteDemandDrivenData(realDataPtr_);
    }
}


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

    interpType_ = interpolationScheme;
    setInterp(interpolationScheme, interpFunc_);
}


template<class Type>
void Foam::lookupTable2D<Type>::setX
(
    const Field<scalar>& x,
    const word& modX,
    const bool isReal
)
{
    modXType_ = modX;
    needXMod_ = modX != "none";

    setMod(modX, modXFunc_, invModXFunc_);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::setX
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
    }
    else
    {
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
    const Field<scalar>& y,
    const word& modY,
    const bool isReal
)
{
    modYType_ = modY;
    needYMod_ = modY != "none";

    setMod(modY, modYFunc_, invModYFunc_);
    setY(y, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::setY
(
    const Field<scalar>& y,
    const bool isReal
)
{
    if (yValuesPtr_ != &yModValues_ && yValuesPtr_ != nullptr)
    {
        deleteDemandDrivenData(yValuesPtr_);
    }

    if (!needYMod_)
    {
        yModValues_ = y;
        yValuesPtr_ = &yModValues_;
    }
    else
    {
        yValuesPtr_ = new scalarField(y);
        yModValues_ = y;

        if (isReal)
        {
            forAll(y, i)
            {
                yModValues_[i] = modYFunc_(y[i]);
            }
        }
        else
        {
            forAll(y, i)
            {
                (*yValuesPtr_)[i] = invModYFunc_(y[i]);
            }
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

    realDataPtr_ = new Field<Field<Type>>(data);
    data_ = data;

    if (isReal)
    {
        forAll(data, i)
        {
            forAll(data[i], j)
            {
                data_[i][j] = modFunc_(data[i][j]);
            }
        }
    }
    else
    {
        forAll(data, i)
        {
            forAll(data[i], j)
            {
                (*realDataPtr_)[i][j] = invModFunc_(data[i][j]);
            }
        }
    }
}


template<class Type>
void Foam::lookupTable2D<Type>::setData
(
    const Field<Field<Type>>& data,
    const word& mod,
    const bool isReal
)
{
    modType_ = mod;
    needMod_ = mod != "none";

    setMod(mod, modFunc_, invModFunc_);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::update(const scalar x, const scalar y) const
{
    scalar xMod(modXFunc_(x));
    scalar yMod(modYFunc_(y));

    ij_.x() = findXIndex_(xMod, xModValues_);
    ij_.y() = findYIndex_(yMod, yModValues_);

    labelList is, js;
    scalarList wxs, wys;

    interpFunc_(xMod, ij_[0], xModValues_, is, wxs);
    interpFunc_(yMod, ij_[1], yModValues_, js, wys);

    indices_.setSize(is.size()*js.size());
    weights_.setSize(indices_.size());

    label n = 0;
    forAll(is, i)
    {
        forAll(js, j)
        {
            indices_[n].x() = is[i];
            indices_[n].y() = js[j];

            weights_[n] = wxs[i]*wys[j];
            n++;
        }
    }
}


template<class Type>
Type Foam::lookupTable2D<Type>::lookup(const scalar x, const scalar y) const
{
    update(x, y);
    Type modf =
        weights_[0]
       *data_
        [indices_[0].x()]
        [indices_[0].y()];

    for (label i = 1; i < indices_.size(); i++)
    {
        modf +=
            weights_[i]
           *data_
            [indices_[i].x()]
            [indices_[i].y()];
    }

    return invModFunc_(modf);
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdX(const scalar x, const scalar y) const
{
    label i = findXIndex_(x, xValues());
    label j = findYIndex_(y, yValues());
    scalar fy = linearWeight(modYFunc_(y), yModValues_[j], yModValues_[j+1]);

    return
        (
            invModFunc_
            (
                data_[i+1][j]*(1.0 - fy)
              + data_[i+1][j+1]*fy
            )
          - invModFunc_
            (
                data_[i][j]*(1.0 - fy)
              + data_[i][j+1]*fy
            )
        )/(xValues()[i+1] - xValues()[i]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdY(const scalar x, const scalar y) const
{
    label i = findXIndex_(x, xValues());
    label j = findYIndex_(y, yValues());
    scalar fx = linearWeight(modXFunc_(x), xModValues_[i], xModValues_[i+1]);

    return
        (
            invModFunc_
            (
                data_[i][j+1]*(1.0 - fx)
              + data_[i+1][j+1]*fx
            )
          - invModFunc_
            (
                data_[i][j]*(1.0 - fx)
              + data_[i+1][j]*fx
            )
        )/(yValues()[j+1] - yValues()[j]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdX2(const scalar x, const scalar y) const
{
    label i = findXIndex_(x, xValues());
    label j = findYIndex_(y, yValues());

    scalar fy = linearWeight(modYFunc_(y), yModValues_[j], yModValues_[j+1]);

    if (i == 0)
    {
        i++;
    }

    const Type gmm(invModFunc_(data_[i-1][j]));
    const Type gm(invModFunc_(data_[i][j]));
    const Type gpm(invModFunc_(data_[i+1][j]));

    const Type gmp(invModFunc_(data_[i-1][j+1]));
    const Type gp(invModFunc_(data_[i][j+1]));
    const Type gpp(invModFunc_(data_[i+1][j+1]));

    const scalar xm(xValues()[i-1]);
    const scalar xi(xValues()[i]);
    const scalar xp(xValues()[i+1]);

    const Type gPrimepm((gpm - gm)/(xp - xi));
    const Type gPrimemm((gm - gmm)/(xi - xm));
    const Type gPrimepp((gpp - gp)/(xp - xi));
    const Type gPrimemp((gp - gmp)/(xi - xm));

    return
        (1.0 - fy)*(gPrimepm - gPrimemm)/(0.5*(xp - xm))
      + fy*(gPrimepp - gPrimemp)/(0.5*(xp - xm));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdY2(const scalar x, const scalar y) const
{
    label i = findXIndex_(x, xValues());
    label j = findYIndex_(y, yValues());
    scalar fx = linearWeight(modXFunc_(x), xModValues_[i], xModValues_[i+1]);

    if (j == 0)
    {
        j++;
    }

    const Type gmm(invModFunc_(data_[i][j-1]));
    const Type gm(invModFunc_(data_[i][j]));
    const Type gmp(invModFunc_(data_[i][j+1]));

    const Type gpm(invModFunc_(data_[i+1][j-1]));
    const Type gp(invModFunc_(data_[i+1][j]));
    const Type gpp(invModFunc_(data_[i+1][j+1]));

    const scalar ym(yValues()[j-1]);
    const scalar yi(yValues()[j]);
    const scalar yp(yValues()[j+1]);

    const Type gPrimemp((gmp - gm)/(yp - yi));
    const Type gPrimemm((gm - gmm)/(yi - ym));
    const Type gPrimepp((gpp - gp)/(yp - yi));
    const Type gPrimepm((gp - gpm)/(yi - ym));

    return
        (1.0 - fx)*(gPrimemp - gPrimemm)/(0.5*(yp - ym))
      + fx*(gPrimepp - gPrimepm)/(0.5*(yp - ym));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdXdY(const scalar x, const scalar y) const
{
    label i = findXIndex_(x, xValues());
    label j = findYIndex_(y, yValues());

    const Type gmm(f()[i][j]);
    const Type gmp(f()[i][j+1]);
    const Type gpm(f()[i+1][j]);
    const Type gpp(f()[i+1][j+1]);

    const scalar xm(xValues()[i]);
    const scalar xp(xValues()[i+1]);

    const scalar ym(yValues()[j]);
    const scalar yp(yValues()[j+1]);

    return ((gpp - gmp)/(xp - xm) - (gpm - gmm)/(xp - xm))/(yp - ym);
}


template<class Type>
void Foam::lookupTable2D<Type>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name,
    const bool canRead
)
{
    interpType_ =
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp");
    setInterp(interpType_, interpFunc_);

    scalarField x;
    bool isReal = readComponent
    (
        dict,
        xName,
        modXType_,
        x
    );
    setX(x, modXType_, isReal);

    scalarField y;
    isReal = readComponent
    (
        dict,
        yName,
        modYType_,
        y
    );
    setY(y, modYType_, isReal);

    const dictionary& fDict(dict.subDict(name + "Coeffs"));
    modType_ = fDict.lookupOrDefault<word>("mod", "none");
    isReal = fDict.lookupOrDefault<Switch>("isReal", true);

    fileName file(fDict.lookup<word>("file"));
    Field<Field<Type>> data
    (
        xModValues_.size(),
        Field<Type>(yModValues_.size())
    );
    read2DTable
    (
        file,
        fDict.lookupOrDefault<string>("delim", ","),
        data,
        fDict.lookupOrDefault<Switch>("flipTable", true),
        !canRead
    );
    setData(data, modType_, isReal);
}

// ************************************************************************* //
