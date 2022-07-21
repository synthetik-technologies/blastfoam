/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "lookupTable3D.H"
#include "tableReader.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

template<class Type>
Type Foam::lookupTable3D<Type>::getValue
(
    const label ijk,
    const scalar f,
    const List<Type>& xyz
) const
{
    if (ijk >= xyz.size())
    {
        return xyz.last();
    }

    return xyz[ijk] + f*(xyz[ijk+1] - xyz[ijk]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable3D<Type>::lookupTable3D()
:
    mod_(nullptr),
    modX_(nullptr),
    modY_(nullptr),
    modZ_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    zInterpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    zValuesPtr_(nullptr),
    ijk_(0, 0, 0),
    indices_(0),
    weights_(0)
{}


template<class Type>
Foam::lookupTable3D<Type>::lookupTable3D(const lookupTable3D<Type>& table)
:
    mod_(table.mod_->clone()),
    modX_(table.modX_->clone()),
    modY_(table.modY_->clone()),
    modZ_(table.modZ_->clone()),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(table.xInterpolator_->clone(xModValues_)),
    yInterpolator_(table.yInterpolator_->clone(yModValues_)),
    zInterpolator_(table.zInterpolator_->clone(zModValues_)),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    zValuesPtr_(nullptr),
    ijk_(0, 0, 0),
    indices_(0),
    weights_(0.0)
{
    set
    (
        table.xModValues_,
        table.yModValues_,
        table.zModValues_,
        table.data_,
        false
    );
}


template<class Type>
Foam::lookupTable3D<Type>::lookupTable3D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& zName,
    const word& name,
    const bool canRead
)
:
    mod_(nullptr),
    modX_(nullptr),
    modY_(nullptr),
    modZ_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    zInterpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    zValuesPtr_(nullptr),
    ijk_(0, 0, 0),
    indices_(0),
    weights_(0)
{
    read(dict, xName, yName, zName, name, canRead);
}


template<class Type>
template<template<class> class ListType1, template<class> class ListType2>
Foam::lookupTable3D<Type>::lookupTable3D
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& z,
    const List<ListType1<ListType2<Type>>>& data,
    const word& modXType,
    const word& modYType,
    const word& modZType,
    const word& modType,
    const word& xInterpolationScheme,
    const word& yInterpolationScheme,
    const word& zInterpolationScheme,
    const bool isReal
)
:
    mod_(Modifier<Type>::New(modType)),
    modX_(Modifier<scalar>::New(modXType)),
    modY_(Modifier<scalar>::New(modYType)),
    modZ_(Modifier<scalar>::New(modZType)),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    zInterpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    yValuesPtr_(nullptr),
    zValuesPtr_(nullptr),
    ijk_(0, 0, 0),
    indices_(0),
    weights_(0)
{
    set
    (
        x, y, z,
        data,
        modXType, modYType, modZType,
        modType,
        xInterpolationScheme,
        yInterpolationScheme,
        zInterpolationScheme,
        isReal
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable3D<Type>::~lookupTable3D()
{
    if (xValuesPtr_ != &xModValues_)
    {
        deleteDemandDrivenData(xValuesPtr_);
    }
    if (yValuesPtr_ != &yModValues_)
    {
        deleteDemandDrivenData(yValuesPtr_);
    }
    if (zValuesPtr_ != &zModValues_)
    {
        deleteDemandDrivenData(zValuesPtr_);
    }
    if (realDataPtr_ != &data_)
    {
        deleteDemandDrivenData(realDataPtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<template<class> class ListType1, template<class> class ListType2>
void Foam::lookupTable3D<Type>::set
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& z,
    const List<ListType1<ListType2<Type>>>& data,
    const bool isReal
)
{
    setX(x, isReal);
    setY(y, isReal);
    setZ(z, isReal);
    setData(data, isReal);
}


template<class Type>
template<template<class> class ListType1, template<class> class ListType2>
void Foam::lookupTable3D<Type>::set
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& z,
    const List<ListType1<ListType2<Type>>>& data,
    const word& modXType,
    const word& modYType,
    const word& modZType,
    const word& modType,
    const word& xInterpolationScheme,
    const word& yInterpolationScheme,
    const word& zInterpolationScheme,
    const bool isReal
)
{
    setX(x, modXType, isReal);
    setY(y, modYType, isReal);
    setZ(z, modZType, isReal);
    setData(data, modType, isReal);

    xInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, xModValues_);
    xInterpolator_->validate();

    yInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, yModValues_);
    yInterpolator_->validate();

    zInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, zModValues_);
    zInterpolator_->validate();
}


template<class Type>
void Foam::lookupTable3D<Type>::setX
(
    const List<scalar>& x,
    const word& modX,
    const bool isReal
)
{
    modX_ = Modifier<scalar>::New(modX);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::setX
(
    const List<scalar>& x,
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
void Foam::lookupTable3D<Type>::setY
(
    const List<scalar>& y,
    const word& modY,
    const bool isReal
)
{
    modY_ = Modifier<scalar>::New(modY);
    setY(y, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::setY
(
    const List<scalar>& y,
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
void Foam::lookupTable3D<Type>::setZ
(
    const List<scalar>& z,
    const word& modZ,
    const bool isReal
)
{
    modZ_ = Modifier<scalar>::New(modZ);
    setZ(z, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::setZ
(
    const List<scalar>& z,
    const bool isReal
)
{
    if (zValuesPtr_ != &zModValues_ && zValuesPtr_ != nullptr)
    {
        deleteDemandDrivenData(zValuesPtr_);
    }

    if (!modZ_->needMod())
    {
        zModValues_ = z;
        zValuesPtr_ = &zModValues_;
    }
    else
    {
        zValuesPtr_ = new scalarField(z);
        zModValues_ = z;

        if (isReal)
        {
            forAll(z, i)
            {
                zModValues_[i] = modZ_()(z[i]);
            }
        }
        else
        {
            forAll(z, i)
            {
                (*zValuesPtr_)[i] = modZ_->inv(z[i]);
            }
        }
    }
    zIndexing_ = indexer::New(zModValues_);
    if (zInterpolator_.valid())
    {
        zInterpolator_->validate();
    }
}


template<class Type>
template<template<class> class ListType1, template<class> class ListType2>
void Foam::lookupTable3D<Type>::setData
(
    const List<ListType1<ListType2<Type>>>& data,
    const bool isReal
)
{
    if (realDataPtr_ != &data_ && realDataPtr_ != nullptr)
    {
        deleteDemandDrivenData(realDataPtr_);
    }

    data_.setSize(data.size());
    forAll(data, i)
    {
        data_[i].setSize(data[i].size());
        forAll(data[i], j)
        {
            data_[i][j] = data[i][j];
        }
    }

    if (!mod_->needMod())
    {
        realDataPtr_ = &data_;
        return;
    }

    realDataPtr_ = new Field<Field<Field<Type>>>(data_);

    if (isReal)
    {
        forAll(data, i)
        {
            forAll(data[i], j)
            {
                forAll(data[i][j], k)
                {
                    data_[i][j][k] = mod_()(data[i][j][k]);
                }
            }
        }
    }
    else
    {
        forAll(data, i)
        {
            forAll(data[i], j)
            {
                forAll(data[i][j], k)
                {
                    (*realDataPtr_)[i][j][k] = mod_->inv(data[i][j][k]);
                }
            }
        }
    }
}


template<class Type>
template<template<class> class ListType1, template<class> class ListType2>
void Foam::lookupTable3D<Type>::setData
(
    const List<ListType1<ListType2<Type>>>& data,
    const word& mod,
    const bool isReal
)
{
    mod_ = Modifier<Type>::New(mod);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::updateIndex
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    ijk_.x() = xIndexing_->findIndex(modX_()(x));
    ijk_.y() = yIndexing_->findIndex(modY_()(y));
    ijk_.z() = zIndexing_->findIndex(modZ_()(z));
}


template<class Type>
void Foam::lookupTable3D<Type>::update
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modX_()(x));
    scalar yMod(modY_()(y));
    scalar zMod(modZ_()(z));

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.y() = yIndexing_->findIndex(yMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    labelList is, js, ks;
    scalarList wxs, wys, wzs;

    xInterpolator_->updateWeights(xMod, ijk_.x(), is, wxs);
    yInterpolator_->updateWeights(yMod, ijk_.y(), js, wys);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks, wzs);

    indices_.setSize(is.size()*js.size()*ks.size());
    weights_.setSize(indices_.size());

    label n = 0;
    forAll(is, i)
    {
        forAll(js, j)
        {
            forAll(ks, k)
            {
                indices_[n].x() = is[i];
                indices_[n].y() = js[j];
                indices_[n].z() = ks[k];

                weights_[n] = wxs[i]*wys[j]*wzs[k];
                n++;
            }
        }
    }
}


template<class Type>
Type Foam::lookupTable3D<Type>::lookup
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    update(x, y, z);
    Type modf =
        weights_[0]
       *data_
        [indices_[0].x()]
        [indices_[0].y()]
        [indices_[0].z()];
    for (label i = 1; i < indices_.size(); i++)
    {
        modf +=
            weights_[i]
           *data_
            [indices_[i].x()]
            [indices_[i].y()]
            [indices_[i].z()];
    }
    return mod_->inv(modf);
}


template<class Type>
Type Foam::lookupTable3D<Type>::dFdX
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar yMod(modY_()(y));
    scalar zMod(modZ_()(z));

    ijk_.x() = xIndexing_->findIndex(modX_()(x));
    const label i = ijk_.x();

    ijk_.y() = yIndexing_->findIndex(yMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    labelList js, ks;
    scalarList wys, wzs;
    yInterpolator_->updateWeights(yMod, ijk_.y(), js, wys);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks, wzs);

    Type fm(data_[i][js[0]][ks[0]]*wys[0]*wzs[0]);
    Type fp(data_[i+1][js[0]][ks[0]]*wys[0]*wzs[0]);
    for (label k = 1; k < ks.size(); k++)
    {
        fm += data_[i][js[0]][ks[k]]*wys[0]*wzs[k];
        fp += data_[i+1][js[0]][ks[k]]*wys[0]*wzs[k];
    }
    for (label j = 1; j < js.size(); j++)
    {
        for (label k = 0; k < ks.size(); k++)
        {
            fm += data_[i][js[j]][ks[k]]*wys[j]*wzs[k];
            fp += data_[i+1][js[j]][ks[k]]*wys[j]*wzs[k];
        }
    }
    return (mod_->inv(fp) - mod_->inv(fm))/(xValues()[i+1] - xValues()[i]);
}


template<class Type>
Type Foam::lookupTable3D<Type>::dFdY
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modX_()(x));
    scalar zMod(modZ_()(z));

    ijk_.y() = yIndexing_->findIndex(modY_()(y));
    const label j = ijk_.y();

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    labelList is, ks;
    scalarList wxs, wzs;
    xInterpolator_->updateWeights(xMod, ijk_.x(), is, wxs);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks, wzs);

    Type fm(data_[is[0]][j][ks[0]]*wxs[0]*wzs[0]);
    Type fp(data_[is[0]][j+1][ks[0]]*wxs[0]*wzs[0]);
    for (label k = 1; k < ks.size(); k++)
    {
        fm += data_[is[0]][j][ks[k]]*wxs[0]*wzs[k];
        fp += data_[is[0]][j+1][ks[k]]*wxs[0]*wzs[k];
    }
    for (label i = 1; i < is.size(); i++)
    {
        for (label k = 0; k < ks.size(); k++)
        {
            fm += data_[is[i]][j][ks[k]]*wxs[i]*wzs[k];
            fp += data_[is[i]][j+1][ks[k]]*wxs[i]*wzs[k];
        }
    }
    return (mod_->inv(fp) - mod_->inv(fm))/(yValues()[j+1] - yValues()[j]);
}

template<class Type>
Type Foam::lookupTable3D<Type>::dFdZ
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modX_()(x));
    scalar yMod(modY_()(y));

    ijk_.z() = zIndexing_->findIndex(modZ_()(z));
    const label k = ijk_.z();

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.y() = yIndexing_->findIndex(yMod);

    labelList is, js;
    scalarList wxs, wys;
    xInterpolator_->updateWeights(xMod, ijk_.x(), is, wxs);
    yInterpolator_->updateWeights(yMod, ijk_.y(), js, wys);

    Type fm(data_[is[0]][js[0]][k]*wxs[0]*wys[0]);
    Type fp(data_[is[0]][js[0]][k+1]*wxs[0]*wys[0]);
    for (label j = 1; j < js.size(); j++)
    {
        fm += data_[is[0]][js[j]][k]*wxs[0]*wys[j];
        fp += data_[is[0]][js[j]][k+1]*wxs[0]*wys[j];
    }
    for (label i = 1; i < is.size(); i++)
    {
        for (label j = 0; j < js.size(); j++)
        {
            fm += data_[is[i]][js[j]][k]*wxs[i]*wys[j];
            fp += data_[is[i]][js[j]][k+1]*wxs[i]*wys[j];
        }
    }
    return (mod_->inv(fp) - mod_->inv(fm))/(zValues()[k+1] - zValues()[k]);
}


template<class Type>
Type Foam::lookupTable3D<Type>::d2FdX2
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar yMod(modY_()(y));
    scalar zMod(modZ_()(z));

    ijk_.x() = xIndexing_->findIndex(modX_()(x));
    const label i = ijk_.x();

    ijk_.y() = yIndexing_->findIndex(yMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    labelList js, ks;
    scalarList wys, wzs;
    yInterpolator_->updateWeights(yMod, ijk_.y(), js, wys);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks, wzs);

    Type fm(data_[i-1][js[0]][ks[0]]*wys[0]*wzs[0]);
    Type f(data_[i][js[0]][ks[0]]*wys[0]*wzs[0]);
    Type fp(data_[i+1][js[0]][ks[0]]*wys[0]*wzs[0]);
    for (label k = 1; k < ks.size(); k++)
    {
        fm += data_[i-1][js[0]][ks[k]]*wys[0]*wzs[k];
        f += data_[i][js[0]][ks[k]]*wys[0]*wzs[k];
        fp += data_[i+1][js[0]][ks[k]]*wys[0]*wzs[k];
    }
    for (label j = 1; j < js.size(); j++)
    {
        for (label k = 0; k < ks.size(); k++)
        {
            fm += data_[i-1][js[j]][ks[k]]*wys[j]*wzs[k];
            f += data_[i][js[j]][ks[k]]*wys[j]*wzs[k];
            fp += data_[i+1][js[j]][ks[k]]*wys[j]*wzs[k];
        }
    }
    const scalar dxm(xValues()[i] - xValues()[i-1]);
    const scalar dxp(xValues()[i+1] - xValues()[i]);

    fm = mod_->inv(fm);
    f = mod_->inv(f);
    fp = mod_->inv(fp);
    return ((fp - f)/dxp - (f - fm)/dxm)/(0.5*(dxp + dxm));
}


template<class Type>
Type Foam::lookupTable3D<Type>::d2FdY2
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modX_()(x));
    scalar zMod(modZ_()(z));

    ijk_.y() = yIndexing_->findIndex(modY_()(y));
    const label j = ijk_.y();

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    labelList is, ks;
    scalarList wxs, wzs;
    xInterpolator_->updateWeights(xMod, ijk_.x(), is, wxs);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks, wzs);

    Type fm(data_[is[0]][j-1][ks[0]]*wxs[0]*wzs[0]);
    Type f(data_[is[0]][j][ks[0]]*wxs[0]*wzs[0]);
    Type fp(data_[is[0]][j+1][ks[0]]*wxs[0]*wzs[0]);
    for (label k = 1; k < ks.size(); k++)
    {
        fm += data_[is[0]][j-1][ks[k]]*wxs[0]*wzs[k];
        f += data_[is[0]][j][ks[k]]*wxs[0]*wzs[k];
        fp += data_[is[0]][j+1][ks[k]]*wxs[0]*wzs[k];
    }
    for (label i = 1; i < is.size(); i++)
    {
        for (label k = 0; k < ks.size(); k++)
        {
            fm += data_[is[i]][j-1][ks[k]]*wxs[i]*wzs[k];
            f += data_[is[i]][j][ks[k]]*wxs[i]*wzs[k];
            fp += data_[is[i]][j+1][ks[k]]*wxs[i]*wzs[k];
        }
    }
    const scalar dym(yValues()[j] - yValues()[j-1]);
    const scalar dyp(yValues()[j+1] - yValues()[j]);
    fm = mod_->inv(fm);
    f = mod_->inv(f);
    fp = mod_->inv(fp);
    return ((fp - f)/dyp - (f - fm)/dym)/(0.5*(dyp + dym));
}

template<class Type>
Type Foam::lookupTable3D<Type>::d2FdZ2
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modX_()(x));
    scalar yMod(modY_()(y));

    ijk_.z() = zIndexing_->findIndex(modZ_()(z));
    const label k = ijk_.z();

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.y() = yIndexing_->findIndex(yMod);

    labelList is, js;
    scalarList wxs, wys;
    xInterpolator_->updateWeights(xMod, ijk_.x(), is, wxs);
    yInterpolator_->updateWeights(yMod, ijk_.y(), js, wys);

    Type fm(data_[is[0]][js[0]][k-1]*wxs[0]*wys[0]);
    Type f(data_[is[0]][js[0]][k]*wxs[0]*wys[0]);
    Type fp(data_[is[0]][js[0]][k+1]*wxs[0]*wys[0]);
    for (label j = 1; j < js.size(); j++)
    {
        fm += data_[is[0]][js[j]][k-1]*wxs[0]*wys[j];
        f += data_[is[0]][js[j]][k]*wxs[0]*wys[j];
        fp += data_[is[0]][js[j]][k+1]*wxs[0]*wys[j];
    }
    for (label i = 1; i < is.size(); i++)
    {
        for (label j = 0; j < js.size(); j++)
        {
            fm += data_[is[i]][js[j]][k-1]*wxs[i]*wys[j];
            f += data_[is[i]][js[j]][k]*wxs[i]*wys[j];
            fp += data_[is[i]][js[j]][k+1]*wxs[i]*wys[j];
        }
    }
    const scalar dzm(zValues()[k] - zValues()[k-1]);
    const scalar dzp(zValues()[k+1] - zValues()[k]);
    fm = mod_->inv(fm);
    f = mod_->inv(f);
    fp = mod_->inv(fp);
    return ((fp - f)/dzp - (f - fm)/dzm)/(0.5*(dzp + dzm));
}


template<class Type>
Type Foam::lookupTable3D<Type>::d2FdXdY
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar zMod(modZ_()(z));

    ijk_.x() = xIndexing_->findIndex(modX_()(x));
    ijk_.y() = yIndexing_->findIndex(modY_()(y));
    const label i = max(ijk_.x(), 1);
    const label j = max(ijk_.y(), 1);

    ijk_.z() = zIndexing_->findIndex(zMod);

    labelList ks;
    scalarList ws;
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks, ws);

    Type fmm(data_[i][j][ks[0]]*ws[0]);
    Type fmp(data_[i][j+1][ks[0]]*ws[0]);
    Type fpm(data_[i+1][j][ks[0]]*ws[0]);
    Type fpp(data_[i+1][j+1][ks[0]]*ws[0]);

    for (label k = 1; k < ks.size(); k++)
    {
        fmm += data_[i][j][ks[k]]*ws[k];
        fmp += data_[i][j+1][ks[k]]*ws[k];
        fpm += data_[i+1][j][ks[k]]*ws[k];
        fpp += data_[i+1][j+1][ks[k]]*ws[k];
    }

    const scalar dx(xValues()[i+1] - xValues()[i]);
    const scalar dy(yValues()[j+1] - yValues()[j]);

    fmm = mod_->inv(fmm);
    fmp = mod_->inv(fmp);
    fpm = mod_->inv(fpm);
    fpp = mod_->inv(fpp);

    return (fpp - fmp - fpm + fmm)/(dx*dy);
}


template<class Type>
Type Foam::lookupTable3D<Type>::d2FdXdZ
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar yMod(modY_()(y));

    ijk_.x() = xIndexing_->findIndex(modX_()(x));
    ijk_.z() = zIndexing_->findIndex(modZ_()(z));
    const label i = max(ijk_.x(), 1);
    const label k = max(ijk_.z(), 1);

    ijk_.y() = yIndexing_->findIndex(yMod);

    labelList js;
    scalarList ws;
    yInterpolator_->updateWeights(yMod, ijk_.y(), js, ws);

    Type fmm(data_[i][js[0]][k]*ws[0]);
    Type fmp(data_[i][js[0]][k+1]*ws[0]);
    Type fpm(data_[i+1][js[0]][k]*ws[0]);
    Type fpp(data_[i+1][js[0]][k+1]*ws[0]);

    for (label j = 1; j < js.size(); j++)
    {
        fmm += data_[i][js[j]][k]*ws[j];
        fmp += data_[i][js[j]][k+1]*ws[j];
        fpm += data_[i+1][js[j]][k]*ws[j];
        fpp += data_[i+1][js[j]][k+1]*ws[j];
    }

    const scalar dx(xValues()[i+1] - xValues()[i]);
    const scalar dz(zValues()[k+1] - zValues()[k]);

    fmm = mod_->inv(fmm);
    fmp = mod_->inv(fmp);
    fpm = mod_->inv(fpm);
    fpp = mod_->inv(fpp);

    return (fpp - fmp - fpm + fmm)/(dx*dz);
}


template<class Type>
Type Foam::lookupTable3D<Type>::d2FdYdZ
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modX_()(x));

    ijk_.y() = yIndexing_->findIndex(modY_()(y));
    ijk_.z() = zIndexing_->findIndex(modZ_()(z));
    const label j = max(ijk_.y(), 1);
    const label k = max(ijk_.z(), 1);

    ijk_.x() = xIndexing_->findIndex(xMod);

    labelList is;
    scalarList ws;
    xInterpolator_->updateWeights(xMod, ijk_.x(), is, ws);

    Type fmm(data_[is[0]][j][k]*ws[0]);
    Type fmp(data_[is[0]][j][k+1]*ws[0]);
    Type fpm(data_[is[0]][j+1][k]*ws[0]);
    Type fpp(data_[is[0]][j+1][k+1]*ws[0]);

    for (label i = 1; i < is.size(); i++)
    {
        fmm += data_[is[i]][j][k]*ws[i];
        fmp += data_[is[i]][j][k+1]*ws[i];
        fpm += data_[is[i]][j+1][k]*ws[i];
        fpp += data_[is[i]][j+1][k+1]*ws[i];
    }

    const scalar dy(yValues()[j+1] - yValues()[j]);
    const scalar dz(zValues()[k+1] - zValues()[k]);

    fmm = mod_->inv(fmm);
    fmp = mod_->inv(fmp);
    fpm = mod_->inv(fpm);
    fpp = mod_->inv(fpp);

    return (fpp - fmp - fpm + fmm)/(dy*dz);
}


template<class Type>
void Foam::lookupTable3D<Type>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& zName,
    const word& name,
    const bool canRead
)
{
    const word scheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );

    scalarField x, y, z;
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

    word modZType;
    isReal = readComponent
    (
        dict,
        zName,
        modZType,
        z
    );
    setZ(z, modZType, isReal);
    zInterpolator_ = interpolationWeight1D::New
    (
        dict.lookupOrDefault<word>
        (
            zName + "InterpolationScheme",
            scheme
        ),
        zModValues_
    );
    zInterpolator_->validate();

    List<List<List<Type>>> data
    (
        xModValues_.size(),
        List<List<Type>>
        (
            yModValues_.size(),
            List<Type>(zModValues_.size())
        )
    );

    word modType = "none";
    if (dict.found(name))
    {
        dict.readIfPresent(name, data);
        dict.readIfPresent(name + "Mod", modType);
        if (modType != "none")
        {
            isReal = dict.lookup<bool>("isReal");
        }
    }
    else if (dict.isDict(name + "Coeffs"))
    {
        const dictionary& fDict(dict.subDict(name + "Coeffs"));
        fDict.readIfPresent("mod", modType);
        if (modType != "none")
        {
            isReal = fDict.lookup<bool>("isReal");
        }

        if (fDict.found(name))
        {
            fDict.readIfPresent(name, data);
        }
        else if (fDict.found("file"))
        {
            fileName file(fDict.lookup<word>("file"));

            read3DTable
            (
                file,
                dict.lookupOrDefault<string>("delim", ","),
                dict.lookupOrDefault<string>("rowDelim", ";"),
                data,
                dict.lookupOrDefault<Switch>("flipTable", true),
                !canRead
            );
        }
        else
        {
            FatalIOErrorInFunction(fDict)
                << "Neither the entry \"" << name << "\", "
                << " or a file was provided for construction" << endl
                << abort(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Neither the entry \"" << name << "\", "
            << " or the \""
            << name << "Coeffs\" subDictionary was found" << endl
            << abort(FatalIOError);
    }

    if
    (
        data.size() != xModValues_.size()
     || data[0].size() != yModValues_.size()
     || data[0][0].size() != zModValues_.size()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Incompatible dimensions for table" << nl
            << "table size: "
            << data.size() << " x "
            << data[0].size() << " x "
            << data[0][0].size() << nl
            << "x and y size: "
            << xModValues_.size() << " x "
            << yModValues_.size() << " z "
            << yModValues_.size() << nl
            << abort(FatalIOError);
    }

    setData(data, modType, isReal);
}

// ************************************************************************* //
