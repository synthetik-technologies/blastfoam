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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable3D<Type>::lookupTable3D()
:
    xName_("x"),
    yName_("y"),
    zName_("z"),
    fName_("f"),
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
    xName_(table.xName_),
    yName_(table.yName_),
    zName_(table.zName_),
    fName_(table.fName_),
    mod_(table.mod_->clone()),
    modX_(table.modX_->clone()),
    modY_(table.modY_->clone()),
    modZ_(table.modZ_->clone()),
    data_(),
    xModValues_(table.xModValues_),
    yModValues_(table.yModValues_),
    zModValues_(table.zModValues_),
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
    xName_(xName),
    yName_(yName),
    zName_(zName),
    fName_(name),
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
template<class ListListListType>
Foam::lookupTable3D<Type>::lookupTable3D
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& z,
    const ListListListType& data,
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
    xName_("x"),
    yName_("y"),
    zName_("z"),
    fName_("f"),
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
template<class ListListListType>
void Foam::lookupTable3D<Type>::set
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& z,
    const ListListListType& data,
    const bool isReal
)
{
    setX(x, isReal);
    setY(y, isReal);
    setZ(z, isReal);
    setData(data, isReal);
}


template<class Type>
template<class ListListListType>
void Foam::lookupTable3D<Type>::set
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& z,
    const ListListListType& data,
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
    xInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, xModValues_);
    xInterpolator_->validate();

    yInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, yModValues_);
    yInterpolator_->validate();

    zInterpolator_ =
        interpolationWeight1D::New(xInterpolationScheme, zModValues_);
    zInterpolator_->validate();

    setX(x, modXType, isReal);
    setY(y, modYType, isReal);
    setZ(z, modZType, isReal);
    setData(data, modType, isReal);
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
        xValuesPtr_ = new scalarList(x);
        xModValues_ = x;

        if (isReal)
        {
            modX_->Mod(xModValues_);
        }
        else
        {
            modX_->Inv(*xValuesPtr_);
        }
    }
    xIndexing_ = indexer::New(xModValues_);
    if (xInterpolator_.valid())
    {
        xInterpolator_->validate();
        xInterpolator_->update();
    }
    else
    {
        xInterpolator_ =
            interpolationWeight1D::New("linearClamp", xModValues_);
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
        yValuesPtr_ = new scalarList(y);
        yModValues_ = y;

        if (isReal)
        {
            modY_->Mod(yModValues_);
        }
        else
        {
            modY_->Inv(*yValuesPtr_);
        }
    }
    yIndexing_ = indexer::New(yModValues_);
    if (yInterpolator_.valid())
    {
        yInterpolator_->validate();
        yInterpolator_->update();
    }
    else
    {
        yInterpolator_ =
            interpolationWeight1D::New("linearClamp", yModValues_);
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
        zValuesPtr_ = new scalarList(z);
        zModValues_ = z;

        if (isReal)
        {
            modZ_->Mod(zModValues_);
        }
        else
        {
            modZ_->Inv(*zValuesPtr_);
        }
    }
    zIndexing_ = indexer::New(zModValues_);
    if (zInterpolator_.valid())
    {
        zInterpolator_->validate();
        zInterpolator_->update();
    }
    else
    {
        zInterpolator_ =
            interpolationWeight1D::New("linearClamp", zModValues_);
        zInterpolator_->validate();
    }
}


template<class Type>
template<class ListListListType>
void Foam::lookupTable3D<Type>::setData
(
    const ListListListType& data,
    const bool isReal
)
{
    if (realDataPtr_ != &data_ && realDataPtr_ != nullptr)
    {
        deleteDemandDrivenData(realDataPtr_);
    }

    data_ = data;

    if (!mod_->needMod())
    {
        realDataPtr_ = &data_;
        return;
    }

    realDataPtr_ = new List3D<Type>(data_);

    if (isReal)
    {
        for (label i = 0; i < data_.m(); i++)
        {
            for (label j = 0; j < data_.n(); j++)
            {
                for (label k = 0; k < data_.l(); k++)
                {
                    data_(i, j, k) = mod()(data[i][j][k]);
                }
            }
        }
    }
    else
    {
        for (label i = 0; i < data_.m(); i++)
        {
            for (label j = 0; j < data_.n(); j++)
            {
                for (label k = 0; k < data_.l(); k++)
                {
                    (*realDataPtr_)(i, j, k) = mod_->inv(data[i][j][k]);
                }
            }
        }
    }
}


template<class Type>
template<class ListListListType>
void Foam::lookupTable3D<Type>::setData
(
    const ListListListType& data,
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

    xInterpolator_->updateWeights(xMod, ijk_.x(), is_, wxs_);
    yInterpolator_->updateWeights(yMod, ijk_.y(), js_, wys_);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks_, wzs_);

    indices_.setSize(is_.size()*js_.size()*ks_.size());
    weights_.setSize(indices_.size());

    label n = 0;
    forAll(is_, i)
    {
        forAll(js_, j)
        {
            forAll(ks_, k)
            {
                indices_[n].x() = is_[i];
                indices_[n].y() = js_[j];
                indices_[n].z() = ks_[k];

                weights_[n] = wxs_[i]*wys_[j]*wzs_[k];
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
    Type modf = weights_[0]*data_(indices_[0]);
    for (label i = 1; i < indices_.size(); i++)
    {
        modf += weights_[i]*data_(indices_[i]);
    }
    return mod_->inv(modf);
}


template<class Type>
Foam::scalar Foam::lookupTable3D<Type>::reverseLookupX
(
    const Type& fin,
    const scalar y,
    const scalar z
) const
{
    NotImplemented;
    return z;
}


template<class Type>
Foam::scalar Foam::lookupTable3D<Type>::reverseLookupY
(
    const Type& fin,
    const scalar x,
    const scalar z
) const
{
    NotImplemented;
    return z;
}


template<class Type>
Foam::scalar Foam::lookupTable3D<Type>::reverseLookupZ
(
    const Type& fin,
    const scalar x,
    const scalar y
) const
{
    NotImplemented;
    return z;
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

    yInterpolator_->updateWeights(yMod, ijk_.y(), js_, wys_);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks_, wzs_);

    Type fm(data_(i, js_[0], ks_[0])*wys_[0]*wzs_[0]);
    Type fp(data_(i+1, js_[0], ks_[0])*wys_[0]*wzs_[0]);
    for (label k = 1; k < ks_.size(); k++)
    {
        fm += data_(i, js_[0], ks_[k])*wys_[0]*wzs_[k];
        fp += data_(i+1, js_[0], ks_[k])*wys_[0]*wzs_[k];
    }
    for (label j = 1; j < js_.size(); j++)
    {
        for (label k = 0; k < ks_.size(); k++)
        {
            fm += data_(i, js_[j], ks_[k])*wys_[j]*wzs_[k];
            fp += data_(i+1, js_[j], ks_[k])*wys_[j]*wzs_[k];
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

    xInterpolator_->updateWeights(xMod, ijk_.x(), is_, wxs_);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks_, wzs_);

    Type fm(data_(is_[0], j, ks_[0])*wxs_[0]*wzs_[0]);
    Type fp(data_(is_[0], j+1, ks_[0])*wxs_[0]*wzs_[0]);
    for (label k = 1; k < ks_.size(); k++)
    {
        fm += data_(is_[0], j, ks_[k])*wxs_[0]*wzs_[k];
        fp += data_(is_[0], j+1, ks_[k])*wxs_[0]*wzs_[k];
    }
    for (label i = 1; i < is_.size(); i++)
    {
        for (label k = 0; k < ks_.size(); k++)
        {
            fm += data_(is_[i], j, ks_[k])*wxs_[i]*wzs_[k];
            fp += data_(is_[i], j+1, ks_[k])*wxs_[i]*wzs_[k];
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

    xInterpolator_->updateWeights(xMod, ijk_.x(), is_, wxs_);
    yInterpolator_->updateWeights(yMod, ijk_.y(), js_, wys_);

    Type fm(data_(is_[0], js_[0], k)*wxs_[0]*wys_[0]);
    Type fp(data_(is_[0], js_[0], k+1)*wxs_[0]*wys_[0]);
    for (label j = 1; j < js_.size(); j++)
    {
        fm += data_(is_[0], js_[j], k)*wxs_[0]*wys_[j];
        fp += data_(is_[0], js_[j], k+1)*wxs_[0]*wys_[j];
    }
    for (label i = 1; i < is_.size(); i++)
    {
        for (label j = 0; j < js_.size(); j++)
        {
            fm += data_(is_[i], js_[j], k)*wxs_[i]*wys_[j];
            fp += data_(is_[i], js_[j], k+1)*wxs_[i]*wys_[j];
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
    const label i = max(ijk_.x(), 1);

    ijk_.y() = yIndexing_->findIndex(yMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    yInterpolator_->updateWeights(yMod, ijk_.y(), js_, wys_);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks_, wzs_);

    Type fm(data_(i-1, js_[0], ks_[0])*wys_[0]*wzs_[0]);
    Type f(data_(i, js_[0], ks_[0])*wys_[0]*wzs_[0]);
    Type fp(data_(i+1, js_[0], ks_[0])*wys_[0]*wzs_[0]);
    for (label k = 1; k < ks_.size(); k++)
    {
        fm += data_(i-1, js_[0], ks_[k])*wys_[0]*wzs_[k];
        f += data_(i, js_[0], ks_[k])*wys_[0]*wzs_[k];
        fp += data_(i+1, js_[0], ks_[k])*wys_[0]*wzs_[k];
    }
    for (label j = 1; j < js_.size(); j++)
    {
        for (label k = 0; k < ks_.size(); k++)
        {
            fm += data_(i-1, js_[j], ks_[k])*wys_[j]*wzs_[k];
            f += data_(i, js_[j], ks_[k])*wys_[j]*wzs_[k];
            fp += data_(i+1, js_[j], ks_[k])*wys_[j]*wzs_[k];
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
    const label j = max(ijk_.y(), 1);

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.z() = zIndexing_->findIndex(zMod);

    xInterpolator_->updateWeights(xMod, ijk_.x(), is_, wxs_);
    zInterpolator_->updateWeights(zMod, ijk_.z(), ks_, wzs_);

    Type fm(data_(is_[0], j-1, ks_[0])*wxs_[0]*wzs_[0]);
    Type f(data_(is_[0], j, ks_[0])*wxs_[0]*wzs_[0]);
    Type fp(data_(is_[0], j+1, ks_[0])*wxs_[0]*wzs_[0]);
    for (label k = 1; k < ks_.size(); k++)
    {
        fm += data_(is_[0], j-1, ks_[k])*wxs_[0]*wzs_[k];
        f += data_(is_[0], j, ks_[k])*wxs_[0]*wzs_[k];
        fp += data_(is_[0], j+1, ks_[k])*wxs_[0]*wzs_[k];
    }
    for (label i = 1; i < is_.size(); i++)
    {
        for (label k = 0; k < ks_.size(); k++)
        {
            fm += data_(is_[i], j-1, ks_[k])*wxs_[i]*wzs_[k];
            f += data_(is_[i], j, ks_[k])*wxs_[i]*wzs_[k];
            fp += data_(is_[i], j+1, ks_[k])*wxs_[i]*wzs_[k];
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
    const label k = max(ijk_.z(), 1);

    ijk_.x() = xIndexing_->findIndex(xMod);
    ijk_.y() = yIndexing_->findIndex(yMod);

    xInterpolator_->updateWeights(xMod, ijk_.x(), is_, wxs_);
    yInterpolator_->updateWeights(yMod, ijk_.y(), js_, wys_);

    Type fm(data_(is_[0], js_[0], k-1)*wxs_[0]*wys_[0]);
    Type f(data_(is_[0], js_[0], k)*wxs_[0]*wys_[0]);
    Type fp(data_(is_[0], js_[0], k+1)*wxs_[0]*wys_[0]);
    for (label j = 1; j < js_.size(); j++)
    {
        fm += data_(is_[0], js_[j], k-1)*wxs_[0]*wys_[j];
        f += data_(is_[0], js_[j], k)*wxs_[0]*wys_[j];
        fp += data_(is_[0], js_[j], k+1)*wxs_[0]*wys_[j];
    }
    for (label i = 1; i < is_.size(); i++)
    {
        for (label j = 0; j < js_.size(); j++)
        {
            fm += data_(is_[i], js_[j], k-1)*wxs_[i]*wys_[j];
            f += data_(is_[i], js_[j], k)*wxs_[i]*wys_[j];
            fp += data_(is_[i], js_[j], k+1)*wxs_[i]*wys_[j];
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

    zInterpolator_->updateWeights(zMod, ijk_.z(), ks_, wzs_);

    Type fmm(data_(i, j, ks_[0])*wzs_[0]);
    Type fmp(data_(i, j+1, ks_[0])*wzs_[0]);
    Type fpm(data_(i+1, j, ks_[0])*wzs_[0]);
    Type fpp(data_(i+1, j+1, ks_[0])*wzs_[0]);

    for (label k = 1; k < ks_.size(); k++)
    {
        fmm += data_(i, j, ks_[k])*wzs_[k];
        fmp += data_(i, j+1, ks_[k])*wzs_[k];
        fpm += data_(i+1, j, ks_[k])*wzs_[k];
        fpp += data_(i+1, j+1, ks_[k])*wzs_[k];
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

    yInterpolator_->updateWeights(yMod, ijk_.y(), js_, wys_);

    Type fmm(data_(i, js_[0], k)*wys_[0]);
    Type fmp(data_(i, js_[0], k+1)*wys_[0]);
    Type fpm(data_(i+1, js_[0], k)*wys_[0]);
    Type fpp(data_(i+1, js_[0], k+1)*wys_[0]);

    for (label j = 1; j < js_.size(); j++)
    {
        fmm += data_(i, js_[j], k)*wys_[j];
        fmp += data_(i, js_[j], k+1)*wys_[j];
        fpm += data_(i+1,js_[j], k)*wys_[j];
        fpp += data_(i+1,js_[j], k+1)*wys_[j];
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

    xInterpolator_->updateWeights(xMod, ijk_.x(), is_, wxs_);

    Type fmm(data_(is_[0], j, k)*wxs_[0]);
    Type fmp(data_(is_[0], j, k+1)*wxs_[0]);
    Type fpm(data_(is_[0], j+1, k)*wxs_[0]);
    Type fpp(data_(is_[0], j+1, k+1)*wxs_[0]);

    for (label i = 1; i < is_.size(); i++)
    {
        fmm += data_(is_[i], j, k)*wxs_[i];
        fmp += data_(is_[i], j, k+1)*wxs_[i];
        fpm += data_(is_[i], j+1, k)*wxs_[i];
        fpp += data_(is_[i], j+1, k+1)*wxs_[i];
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
    xName_ = xName;
    yName_ = yName;
    zName_ = zName;
    fName_ = name;

    const word scheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );

    scalarField x;
    {
        const dictionary& xDict = readComponent
        (
            dict,
            xName,
            modX_,
            x,
            canRead
        );
        setX(x, true);
        xInterpolator_ = interpolationWeight1D::New
        (
            xDict.found("interpolationScheme")
          ? xDict.lookup<word>("interpolationScheme")
          : dict.lookupOrDefault<word>
            (
                xName + "InterpolationScheme",
                scheme
            ),
            xModValues_,
            canRead
        );
        xInterpolator_->validate();
    }

    scalarField y;
    {
        const dictionary& yDict = readComponent
        (
            dict,
            yName,
            modY_,
            y,
            canRead
        );
        setY(y, true);
        yInterpolator_ = interpolationWeight1D::New
        (
            yDict.found("interpolationScheme")
          ? yDict.lookup<word>("interpolationScheme")
          : dict.lookupOrDefault<word>
            (
                yName + "InterpolationScheme",
                scheme
            ),
            yModValues_,
            canRead
        );
        yInterpolator_->validate();
    }

    scalarField z;
    {
        const dictionary& zDict = readComponent
        (
            dict,
            zName,
            modZ_,
            z,
            canRead
        );
        setZ(z, true);
        zInterpolator_ = interpolationWeight1D::New
        (
            zDict.found("interpolationScheme")
          ? zDict.lookup<word>("interpolationScheme")
          : dict.lookupOrDefault<word>
            (
                zName + "InterpolationScheme",
                scheme
            ),
            zModValues_,
            canRead
        );
        zInterpolator_->validate();
    }

    List3D<Type> data
    (
        xModValues_.size(),
        yModValues_.size(),
        zModValues_.size()
    );

    if (dict.found(name))
    {
        dict.readIfPresent(name, data);
        mod_ = Modifier<Type>::New
        (
            dict.lookup<word>(name + "Mod"),
            dict
        );
        mod_->readReal(dict, name + "IsReal");
    }
    else if (dict.isDict(name + "Coeffs"))
    {
        const dictionary& fDict(dict.subDict(name + "Coeffs"));
        mod_ = Modifier<Type>::New(fDict.lookup<word>("mod"), fDict);
        mod_->readReal(fDict, "isReal");

        if (fDict.found(name))
        {
            fDict.readIfPresent(name, data);
        }
        else if (fDict.found("file"))
        {
            read3DTable
            (
                fDict.lookup<fileName>("file"),
                readDelim(dict),
                readDelim(dict, "rowDelim", token::END_STATEMENT),
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
        data.m() != xModValues_.size()
     || data.n() != yModValues_.size()
     || data.l() != zModValues_.size()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Incompatible dimensions for table" << nl
            << "table size: "
            << data.m() << " x "
            << data.n() << " x "
            << data.l() << nl
            << "x and y size: "
            << xModValues_.size() << " x "
            << yModValues_.size() << " z "
            << yModValues_.size() << nl
            << abort(FatalIOError);
    }
    if (!mod_->isReal())
    {
        mod_->Inv(data);
        mod_->setReal();
    }
    setData(data, true);

    if (dict.found("rootSolver"))
    {
        this->solver
        (
            dict.lookup<word>("rootSolver"),
            dict
        );
    }
}


template<class Type>
void  Foam::lookupTable3D<Type>::write(Ostream& os) const
{
    os << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    if (solver_.valid())
    {
        writeEntry(os, "rootSolver", solver_->type());
    }

    writeEntry(os, word(xName_ + "InterpolationScheme"), xInterpolator_->type());
    writeEntry(os, word(yName_ + "InterpolationScheme"), yInterpolator_->type());
    writeEntry(os, word(zName_ + "InterpolationScheme"), zInterpolator_->type());

    writeKeyword(os, word(xName_ + "Coeffs"))
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

        writeEntry<bool>(os, "isReal", true);
        writeEntry(os, "mod", modX_->type());
        writeEntry(os, xName_, static_cast<const scalarList&>(x()));

    os  << decrIndent << indent << token::END_BLOCK << endl;

    writeKeyword(os, word(yName_ + "Coeffs"))
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

        writeEntry<bool>(os, "isReal", true);
        writeEntry(os, "mod", modY_->type());
        writeEntry(os, yName_, static_cast<const scalarList&>(y()));

    os  << decrIndent << indent << token::END_BLOCK << endl;

    writeKeyword(os, word(zName_ + "Coeffs"))
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

        writeEntry<bool>(os, "isReal", true);
        writeEntry(os, "mod", modZ_->type());
        writeEntry(os, zName_, static_cast<const scalarList&>(z()));

    os  << decrIndent << indent << token::END_BLOCK << endl;


    writeKeyword(os, word(fName_ + "Coeffs"))
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

        writeEntry<bool>(os, "isReal", true);
        writeEntry(os, "mod", mod_->type());
        writeEntry(os, fName_, f());

    os  << decrIndent << indent << token::END_BLOCK << endl;

    os  << decrIndent << indent << token::END_BLOCK << endl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable3D<Type>::operator=(const lookupTable3D<Type>& table)
{
    xName_ = table.xName_;
    yName_ = table.yName_;
    zName_ = table.zName_;
    fName_ = table.fName_;

    mod_ = table.mod_->clone();
    modX_ = table.modX_->clone();
    modY_ = table.modY_->clone();
    modZ_ = table.modZ_->clone();

    xInterpolator_ = table.xInterpolator_->clone(xModValues_);
    yInterpolator_ = table.yInterpolator_->clone(yModValues_);
    zInterpolator_ = table.zInterpolator_->clone(zModValues_);
    set
    (
        table.xModValues_,
        table.yModValues_,
        table.zModValues_,
        table.data_,
        false
    );
}

// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class Type>
void  Foam::writeEntry(Ostream& os, const lookupTable3D<Type>& table)
{
    table.write(os);
}

// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const lookupTable3D<Type>& table
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const lookupTable3D<Type>&)"
    );

    table.write(os);

    return os;
}

// ************************************************************************* //
