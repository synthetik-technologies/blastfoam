/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D()
:
    xName_("x"),
    yName_("y"),
    fName_("f"),
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
    xName_(table.xName_),
    yName_(table.yName_),
    fName_(table.fName_),
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
    xName_(xName),
    yName_(yName),
    fName_(name),
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
template<template<class> class ListListType>
Foam::lookupTable2D<Type>::lookupTable2D
(
    const List<scalar>& x,
    const List<scalar>& y,
    const ListListType<Type>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const word& xInterpolationScheme,
    const word& yInterpolationScheme,
    const bool isReal
)
:
    xName_("x"),
    yName_("y"),
    fName_("f"),
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
template<template<class> class ListListType>
void Foam::lookupTable2D<Type>::set
(
    const List<scalar>& x,
    const List<scalar>& y,
    const ListListType<Type>& data,
    const bool isReal
)
{
    setX(x, isReal);
    setY(y, isReal);
    setData(data, isReal);
}


template<class Type>
template<template<class> class ListListType>
void Foam::lookupTable2D<Type>::set
(
    const List<scalar>& x,
    const List<scalar>& y,
    const ListListType<Type>& data,
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
    const List<scalar>& x,
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
    const List<scalar>& y,
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
template<template<class> class ListListType>
void Foam::lookupTable2D<Type>::setData
(
    const ListListType<Type>& data,
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

    realDataPtr_ = new List2D<Type>(data_);
    if (isReal)
    {
        for (label i = 0; i < data_.m(); i++)
        {
            for (label j = 0; j < data_.n(); j++)
            {
                data_(i, j) = mod_()(data[i][j]);
            }
        }
    }
    else
    {
        for (label i = 0; i < data_.m(); i++)
        {
            for (label j = 0; j < data_.n(); j++)
            {
                (*realDataPtr_)(i, j) = mod_->inv(data[i][j]);
            }
        }
    }
}


template<class Type>
template<template<class> class ListListType>
void Foam::lookupTable2D<Type>::setData
(
    const ListListType<Type>& data,
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

    xInterpolator_->updateWeights(xMod, ij_[0], is_, wxs_);
    yInterpolator_->updateWeights(yMod, ij_[1], js_, wys_);

    indices_.setSize(is_.size()*js_.size());
    weights_.setSize(indices_.size());

    label n = 0;
    forAll(is_, i)
    {
        forAll(js_, j)
        {
            indices_[n].x() = is_[i];
            indices_[n].y() = js_[j];

            weights_[n] = wxs_[i]*wys_[j];
            n++;
        }
    }
}


template<class Type>
Type Foam::lookupTable2D<Type>::lookup(const scalar x, const scalar y) const
{
    update(x, y);
    Type modf = weights_[0]*data_(indices_[0]);
    for (label i = 1; i < indices_.size(); i++)
    {
        modf += weights_[i]*data_(indices_[i]);
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
    NotImplemented
    return y;
}


template<class Type>
Foam::scalar Foam::lookupTable2D<Type>::reverseLookupY
(
    const Type& fin,
    const scalar x
) const
{
    NotImplemented;
    return x;
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdX(const scalar x, const scalar y) const
{
    scalar yMod(modY_()(y));
    ij_.x() = xIndexing_->findIndex(modX_()(x));
    const label i = ij_.x();

    ij_.y() = yIndexing_->findIndex(yMod);

    yInterpolator_->updateWeights(yMod, ij_.y(), js_, wys_);
    Type fm(data_(i, js_[0])*wys_[0]);
    Type fp(data_(i+1, js_[0])*wys_[0]);
    for (label j = 1; j < js_.size(); j++)
    {
        fm += data_(i, js_[j])*wys_[j];
        fp += data_(i+1, js_[j])*wys_[j];
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

    xInterpolator_->updateWeights(xMod, ij_.x(), is_, wxs_);
    Type fm(data_(is_[0], j)*wxs_[0]);
    Type fp(data_(is_[0], j+1)*wxs_[0]);
    for (label i = 1; i < is_.size(); i++)
    {
        fm += data_(is_[i], j)*wxs_[i];
        fp += data_(is_[i], j+1)*wxs_[i];
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

    yInterpolator_->updateWeights(yMod, ij_.y(), js_, wys_);

    const scalar dxm(xValues()[i] - xValues()[i-1]);
    const scalar dxp(xValues()[i+1] - xValues()[i]);

    Type fm(data_(i-1, js_[0])*wys_[0]);
    Type f(data_(i, js_[0])*wys_[0]);
    Type fp(data_(i+1, js_[0])*wys_[0]);
    for (label j = 1; j < js_.size(); j++)
    {
        fm += data_(i-1, js_[j])*wys_[j];
        f += data_(i, js_[j])*wys_[j];
        fp += data_(i+1, js_[j])*wys_[j];
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

    const scalar dym(yValues()[j] - yValues()[j-1]);
    const scalar dyp(yValues()[j+1] - yValues()[j]);

    xInterpolator_->updateWeights(xMod, ij_.x(), is_, wxs_);
    Type fm(data_(is_[0], j-1)*wxs_[0]);
    Type f(data_(is_[0], j)*wxs_[0]);
    Type fp(data_(is_[0], j+1)*wxs_[0]);
    for (label i = 1; i < is_.size(); i++)
    {
        fm += data_(is_[i], j-1)*wxs_[i];
        f += data_(is_[i], j)*wxs_[i];
        fp += data_(is_[i], j+1)*wxs_[i];
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

    const Type& fmm(f()(i, j));
    const Type& fmp(f()(i, j+1));
    const Type& fpm(f()(i+1, j));
    const Type& fpp(f()(i+1, j+1));

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
    xName_ = xName;
    yName_ = yName;
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
        if (canRead)
        {
            xInterpolator_->validate();
        }
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
        if (canRead)
        {
            yInterpolator_->validate();
        }
    }

    List2D<Type> data(xModValues_.size(), yModValues_.size());
    if (dict.found(name) || dict.found(name + "File"))
    {
        if (canRead)
        {
            if (dict.found(name))
            {
                dict.readIfPresent(name, data);
            }
            else
            {
                read2DTable
                (
                    dict.lookup<fileName>(name + "File"),
                    dict.lookupOrDefault<char>(name + "Delim", ','),
                    data,
                    dict.lookupOrDefault<bool>(name + "FlipTable", false),
                    !canRead
                );
            }
        }
        mod_ = Modifier<scalar>::New
        (
            dict.lookupOrDefault<word>(name + "Mod", "none"),
            dict
        );
        mod_->readReal(dict, name + "IsReal");
    }
    else if (dict.isDict(name + "Coeffs"))
    {
        const dictionary& fDict(dict.subDict(name + "Coeffs"));
        mod_ = Modifier<scalar>::New
        (
            fDict.lookupOrDefault<word>("mod", "none"),
            fDict
        );
        mod_->readReal(fDict, "isReal");

        if (!canRead)
        {}
        else if (fDict.found(name))
        {
            fDict.readIfPresent(name, data);
        }
        else if (fDict.found("file"))
        {
            read2DTable
            (
                fDict.lookup<fileName>("file"),
                fDict.lookupOrDefault<char>("delim", ','),
                data,
                fDict.lookupOrDefault<bool>("flipTable", false),
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
    else if (canRead)
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
    )
    {
        FatalIOErrorInFunction(dict)
            << "Incompatible dimensions for table" << nl
            << "table size: "
            << data.m() << " x " << data.n() << nl
            << "x and y size: "
            << xModValues_.size() << " x " << yModValues_.size() << nl
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
void  Foam::lookupTable2D<Type>::write(Ostream& os) const
{
    os << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    if (solver_.valid())
    {
        writeEntry(os, "rootSolver", solver_->type());
    }

    writeEntry(os, word(xName_ + "InterpolationScheme"), xInterpolator_->type());
    writeEntry(os, word(yName_ + "InterpolationScheme"), yInterpolator_->type());

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


    writeKeyword(os, word(fName_ + "Coeffs"))
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

        writeEntry<bool>(os, "isReal", true);
        writeEntry(os, "mod", mod_->type());
        writeEntry(os, fName_, this->f());

    os  << decrIndent << indent << token::END_BLOCK << endl;

    os  << decrIndent << indent << token::END_BLOCK << endl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable2D<Type>::operator=(const lookupTable2D<Type>& table)
{
    xName_ = table.xName_;
    yName_ = table.yName_;
    fName_ = table.fName_;

    mod_ = table.mod_->clone();
    modX_ = table.modX_->clone();
    modY_ = table.modY_->clone();

    xInterpolator_ = table.xInterpolator_->clone(xModValues_);
    yInterpolator_ = table.yInterpolator_->clone(yModValues_);
    set(table.xModValues_, table.yModValues_, table.data_, false);
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class Type>
void  Foam::writeEntry(Ostream& os, const lookupTable2D<Type>& table)
{
    table.write(os);
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const lookupTable2D<Type>& table
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const lookupTable2D<Type>&)"
    );

    table.write(os);

    return os;
}

// ************************************************************************* //
