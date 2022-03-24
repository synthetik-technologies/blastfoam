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
Foam::labelList Foam::lookupTable2D<Type>::boundi(const Type& f) const
{
    NotImplemented;
    return labelList();
}

template<class Type>
Foam::labelList Foam::lookupTable2D<Type>::boundj(const Type& f) const
{
    NotImplemented;
    return labelList();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D()
:
    mod_(nullptr),
    modX_(nullptr),
    modY_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    ij_(0, 0),
    indices_(0),
    weights_(0)
{}


template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D(const lookupTable2D<Type>& table)
:
    mod_(table.mod_->clone()),
    modX_(table.modX_->clone()),
    modY_(table.modY_->clone()),
    data_(),
    xModValues_(),
    yModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    xInterpolator_(table.xInterpolator_->clone(xModValues_)),
    yInterpolator_(table.yInterpolator_->clone(yModValues_)),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    ij_(0, 0),
    indices_(0),
    weights_(0.0)
{
    set
    (
        table.xModValues_,
        table.yModValues_,
        table.data_,
        false
    );
}


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
    mod_(nullptr),
    modX_(nullptr),
    modY_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
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
    const word& xInterpolationScheme,
    const word& yInterpolationScheme,
    const bool isReal
)
:
    mod_(Modifier<Type>::New(modType)),
    modX_(Modifier<scalar>::New(modXType)),
    modY_(Modifier<scalar>::New(modYType)),
    data_(data),
    xModValues_(x),
    yModValues_(y),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    ij_(0, 0),
    indices_(0),
    weights_(0)
{
    set
    (
        x,
        y,
        data,
        modXType,
        modYType,
        modType,
        xInterpolationScheme,
        yInterpolationScheme,
        isReal
    );
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
    const word& xInterpolationScheme,
    const word& yInterpolationScheme,
    const bool isReal
)
{
    setX(x, modXType, isReal);
    setY(y, modYType, isReal);
    setData(data, modType, isReal);

    xInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, xModValues_);
    xInterpolator_->validate();

    yInterpolator_ =
        interpolationWeight1D::New(yInterpolationScheme, yModValues_);
    yInterpolator_->validate();
}


template<class Type>
void Foam::lookupTable2D<Type>::setX
(
    const Field<scalar>& x,
    const word& modX,
    const bool isReal
)
{
    modX_ = Modifier<scalar>::New(modX);
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

    if (!modX_->needMod())
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
                xModValues_[i] = modX_()(x[i]);
            }
        }
        else
        {
            forAll(x, i)
            {
                (*xValuesPtr_)[i] = modX_->inv(x[i]);
            }
        }
    }
    xIndexing_ = indexer::New(xModValues_);
    if (xInterpolator_.valid())
    {
        xInterpolator_->validate();
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
    modY_ = Modifier<scalar>::New(modY);
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

    if (!modY_->needMod())
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
                yModValues_[i] = modY_()(y[i]);
            }
        }
        else
        {
            forAll(y, i)
            {
                (*yValuesPtr_)[i] = modY_->inv(y[i]);
            }
        }
    }
    yIndexing_ = indexer::New(yModValues_);
    if (yInterpolator_.valid())
    {
        yInterpolator_->validate();
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

    if (!mod_->needMod())
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
                data_[i][j] = mod_()(data[i][j]);
            }
        }
    }
    else
    {
        forAll(data, i)
        {
            forAll(data[i], j)
            {
                (*realDataPtr_)[i][j] = mod_->inv(data[i][j]);
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
    mod_ = Modifier<Type>::New(mod);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::updateIndex
(
    const scalar x,
    const scalar y
) const
{
    ij_.x() = xIndexing_->findIndex(modX_()(x));
    ij_.y() = yIndexing_->findIndex(modY_()(y));
}


template<class Type>
void Foam::lookupTable2D<Type>::update(const scalar x, const scalar y) const
{
    scalar xMod(modX_()(x));
    scalar yMod(modY_()(y));

    ij_.x() = xIndexing_->findIndex(xMod);
    ij_.y() = yIndexing_->findIndex(yMod);

    labelList is, js;
    scalarList wxs, wys;

    xInterpolator_->updateWeights(xMod, ij_[0], is, wxs);
    yInterpolator_->updateWeights(yMod, ij_[1], js, wys);

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

    return mod_->inv(modf);
}

template<class Type>
Foam::scalar Foam::lookupTable2D<Type>::reverseLookupX
(
    const Type& fin,
    const scalar y
) const
{
    NotImplemented;
    return 0.0;
}


template<class Type>
Foam::scalar Foam::lookupTable2D<Type>::reverseLookupY
(
    const Type& fin,
    const scalar x
) const
{
    NotImplemented;
    return 0.0;
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdX(const scalar x, const scalar y) const
{
    scalar yMod(modY_()(y));
    ij_.x() = xIndexing_->findIndex(modX_()(x));
    const label i = ij_.x();

    ij_.y() = yIndexing_->findIndex(yMod);

    labelList js;
    scalarList ws;
    yInterpolator_->updateWeights(yMod, ij_.y(), js, ws);
    Type fm(data_[i][js[0]]*ws[0]);
    Type fp(data_[i+1][js[0]]*ws[0]);
    for (label j = 1; j < js.size(); j++)
    {
        fm += data_[i][js[j]]*ws[j];
        fp += data_[i+1][js[j]]*ws[j];
    }

    return
        (mod_->inv(fp) - mod_->inv(fm))/(xValues()[i+1] - xValues()[i]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdY(const scalar x, const scalar y) const
{
    scalar xMod(modX_()(x));
    ij_.x() = xIndexing_->findIndex(xMod);
    ij_.y() = yIndexing_->findIndex(modY_()(y));
    const label j = ij_.y();

    labelList is;
    scalarList ws;
    xInterpolator_->updateWeights(xMod, ij_.x(), is, ws);
    Type fm(data_[is[0]][j]*ws[0]);
    Type fp(data_[is[0]][j+1]*ws[0]);
    for (label i = 1; i < is.size(); i++)
    {
        fm += data_[is[i]][j]*ws[i];
        fp += data_[is[i]][j+1]*ws[i];
    }

    return (mod_->inv(fp) - mod_->inv(fm))/(yValues()[j+1] - yValues()[j]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdX2(const scalar x, const scalar y) const
{
    scalar yMod(modY_()(y));
    ij_.x() = xIndexing_->findIndex(modX_()(x));
    const label i = max(ij_.x(), 1);

    ij_.y() = yIndexing_->findIndex(yMod);

    labelList js;
    scalarList ws;
    yInterpolator_->updateWeights(yMod, ij_.y(), js, ws);

    const scalar dxm(xValues()[i] - xValues()[i-1]);
    const scalar dxp(xValues()[i+1] - xValues()[i]);

    Type fm(data_[i-1][js[0]]*ws[0]);
    Type f(data_[i][js[0]]*ws[0]);
    Type fp(data_[i+1][js[0]]*ws[0]);
    for (label j = 1; j < js.size(); j++)
    {
        fm += data_[i-1][js[j]]*ws[j];
        f += data_[i][js[j]]*ws[j];
        fp += data_[i+1][js[j]]*ws[j];
    }
    fm = mod_->inv(fm);
    f = mod_->inv(f);
    fp = mod_->inv(fp);
    return ((fp - f)/dxp - (f - fm)/dxm)/(0.5*(dxp + dxm));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdY2(const scalar x, const scalar y) const
{
    scalar xMod(modX_()(x));
    ij_.x() = xIndexing_->findIndex(xMod);
    ij_.y() = yIndexing_->findIndex(modY_()(y));
    const label j = max(ij_.y(), 1);

    labelList is;
    scalarList ws;

    const scalar dym(yValues()[j] - yValues()[j-1]);
    const scalar dyp(yValues()[j+1] - yValues()[j]);

    xInterpolator_->updateWeights(xMod, ij_.x(), is, ws);
    Type fm(data_[is[0]][j-1]*ws[0]);
    Type f(data_[is[0]][j]*ws[0]);
    Type fp(data_[is[0]][j+1]*ws[0]);
    for (label i = 1; i < is.size(); i++)
    {
        fm += data_[is[i]][j-1]*ws[i];
        f += data_[is[i]][j]*ws[i];
        fp += data_[is[i]][j+1]*ws[i];
    }

    fm = mod_->inv(fm);
    f = mod_->inv(f);
    fp = mod_->inv(fp);
    return ((fp - f)/dyp - (f - fm)/dym)/(0.5*(dyp + dym));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdXdY(const scalar x, const scalar y) const
{
    label i = xIndexing_->findIndex(modX_()(x));
    label j = yIndexing_->findIndex(modY_()(y));

    const Type& fmm(f()[i][j]);
    const Type& fmp(f()[i][j+1]);
    const Type& fpm(f()[i+1][j]);
    const Type& fpp(f()[i+1][j+1]);

    const scalar xm(xValues()[i]);
    const scalar xp(xValues()[i+1]);

    const scalar ym(yValues()[j]);
    const scalar yp(yValues()[j+1]);

    return ((fpp - fmp)/(xp - xm) - (fpm - fmm)/(xp - xm))/(yp - ym);
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
    const word scheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );

    scalarField x;
    word modXType;
    bool isReal = readComponent
    (
        dict,
        xName,
        modXType,
        x
    );
    setX(x, modXType, isReal);
    xInterpolator_ = interpolationWeight1D::New
    (
        dict.lookupOrDefault<word>
        (
            xName + "InterpolationScheme",
            scheme
        ),
        xModValues_
    );
    xInterpolator_->validate();

    scalarField y;
    word modYType;
    isReal = readComponent
    (
        dict,
        yName,
        modYType,
        y
    );
    setY(y, modYType, isReal);
    yInterpolator_ = interpolationWeight1D::New
    (
        dict.lookupOrDefault<word>
        (
            yName + "InterpolationScheme",
            scheme
        ),
        yModValues_
    );
    yInterpolator_->validate();

    const dictionary& fDict(dict.subDict(name + "Coeffs"));
    word modType = fDict.lookupOrDefault<word>("mod", "none");
    if (modType != "none")
    {
        isReal = fDict.lookup<Switch>("isReal");
    }

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
    setData(data, modType, isReal);
}

// ************************************************************************* //
