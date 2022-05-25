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
    mod_(Modifier<Type>::New("none")),
    modX_(Modifier<scalar>::New("none")),
    data_(),
    xModValues_(),
    indexing_(nullptr),
    interpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0)
{}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D(const lookupTable1D<Type>& table)
:
    mod_(table.mod_->clone()),
    modX_(table.modX_->clone()),
    data_(),
    xModValues_(),
    indexing_(nullptr),
    interpolator_(table.interpolator_->clone(xModValues_)),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0)
{
    set(table.xModValues_, table.data_, false);
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
    mod_(nullptr),
    modX_(nullptr),
    data_(),
    xModValues_(),
    indexing_(nullptr),
    interpolator_(nullptr),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0)
{
    read(dict, xName, name, canRead);
}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D
(
    const List<scalar>& x,
    const List<Type>& data,
    const word& xMod,
    const word& mod,
    const word& interpolationScheme,
    const bool isReal
)
:
    mod_(Modifier<Type>::New(mod)),
    modX_(Modifier<scalar>::New(xMod)),
    data_(data),
    xModValues_(x),
    indexing_(nullptr),
    interpolator_(interpolationWeight1D::New(interpolationScheme, xModValues_)),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0)
{
    set(x, data, isReal);
}


template<class Type>
Foam::lookupTable1D<Type>::lookupTable1D
(
    const List<scalar>& x,
    const word& xMod,
    const word& interpolationScheme,
    const bool isReal
)
:
    mod_(Modifier<Type>::New("none")),
    modX_(Modifier<scalar>::New(xMod)),
    data_(),
    xModValues_(),
    indexing_(nullptr),
    interpolator_(interpolationWeight1D::New(interpolationScheme, xModValues_)),
    realDataPtr_(nullptr),
    xValuesPtr_(nullptr),
    index_(0),
    indices_(0),
    weights_(0)
{
    setX(x, isReal);
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
    const List<scalar>& x,
    const List<Type>& data,
    const bool isReal
)
{
    setX(x, isReal);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::set
(
    const List<scalar>& x,
    const List<Type>& data,
    const word& xMod,
    const word& mod,
    const word& interpolationScheme,
    const bool isReal
)
{
    setX(x, xMod, isReal);
    setData(data, mod, isReal);
    interpolator_ =
        interpolationWeight1D::New(interpolationScheme, xModValues_);

}


template<class Type>
void Foam::lookupTable1D<Type>::setX
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
                (*xValuesPtr_)[i] = modX_().inv(x[i]);
            }
        }
    }
    indexing_ = indexer::New(xModValues_);
    if (interpolator_.valid())
    {
        interpolator_->validate();
    }
}


template<class Type>
void Foam::lookupTable1D<Type>::setX
(
    const List<scalar>& x,
    const word& xMod,
    const bool isReal
)
{
    modX_ = Modifier<scalar>::New(xMod);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::setData
(
    const List<Type>& data,
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

    realDataPtr_ = new Field<Type>(data);
    data_ = data;
    if (isReal)
    {
        forAll(data, i)
        {
            data_[i] = mod_()(data[i]);
        }
    }
    else
    {
        forAll(data, i)
        {
            (*realDataPtr_)[i] = mod_->inv(data[i]);
        }
    }
}


template<class Type>
void Foam::lookupTable1D<Type>::setData
(
    const List<Type>& data,
    const word& mod,
    const bool isReal
)
{
    mod_ = Modifier<Type>::New(mod);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable1D<Type>::updateIndex(const scalar x) const
{
    index_ = indexing_->findIndex(modX_()(x));
}


template<class Type>
void Foam::lookupTable1D<Type>::update(const scalar x) const
{
    scalar xMod(modX_()(x));
    index_ = indexing_->findIndex(xMod);
    interpolator_->updateWeights
    (
        xMod,
        index_,
        indices_,
        weights_
    );
}


template<class Type>
Type Foam::lookupTable1D<Type>::lookup(const scalar x) const
{
#ifdef FULL_DEBUG
    if (!mod_.valid())
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
    return mod_->inv(modf);
}


template<class Type>
Type Foam::lookupTable1D<Type>::dFdX(const scalar x) const
{
#ifdef FULL_DEBUG
    if (!mod_.valid())
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    updateIndex(x);

    return
        (f()[index_ + 1] - f()[index_])
       /(xValues()[index_ + 1] - xValues()[index_]);
}


template<class Type>
Type Foam::lookupTable1D<Type>::d2FdX2(const scalar x) const
{
#ifdef FULL_DEBUG
    if (!mod_.valid())
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    updateIndex(x);
    if (index_ == 0)
    {
        index_++;
    }

    const Type& ym(f()[index_-1]);
    const Type& yi(f()[index_]);
    const Type& yp(f()[index_+1]);

    const scalar dxm(xValues()[index_] - xValues()[index_-1]);
    const scalar dxp(xValues()[index_+1] - xValues()[index_]);

    return ((yp - yi)/dxp - (yi - ym)/dxm)/(0.5*(dxp + dxm));
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
    word modXType;
    bool isReal = readComponent<scalar>
    (
        dict,
        xName,
        modXType,
        x,
        table
    );
    setX(x, modXType, isReal);
    interpolator_ =
        interpolationWeight1D::New
        (
            dict.lookupOrDefault<word>("interpolationScheme", "linearClamp"),
            xModValues_
        );
    interpolator_->validate();

    Field<Type> data;
    word modType;
    isReal = readComponent<Type>
    (
        dict,
        name,
        modType,
        data,
        table
    );
    setData(data, modType, isReal);
}

// ************************************************************************* //
