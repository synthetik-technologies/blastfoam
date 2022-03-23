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
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    zInterpolator_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
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
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(table.xInterpolator_->clone()),
    yInterpolator_(table.yInterpolator_->clone()),
    zInterpolator_(table.zInterpolator_->clone()),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
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
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    zInterpolator_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
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
Foam::lookupTable3D<Type>::lookupTable3D
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<scalar>& z,
    const Field<Field<Field<Type>>>& data,
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
    xIndexing_(nullptr),
    yIndexing_(nullptr),
    zIndexing_(nullptr),
    xInterpolator_(nullptr),
    yInterpolator_(nullptr),
    zInterpolator_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
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
void Foam::lookupTable3D<Type>::set
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<scalar>& z,
    const Field<Field<Field<Type>>>& data,
    const bool isReal
)
{
    setX(x, isReal);
    setY(y, isReal);
    setZ(z, isReal);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::set
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<scalar>& z,
    const Field<Field<Field<Type>>>& data,
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

    xInterpolator_ = interpolationWeight1D::New(xInterpolationScheme, x.size());
    yInterpolator_ = interpolationWeight1D::New(xInterpolationScheme, y.size());
    zInterpolator_ = interpolationWeight1D::New(xInterpolationScheme, z.size());
}


template<class Type>
void Foam::lookupTable3D<Type>::setX
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
void Foam::lookupTable3D<Type>::setX
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
}


template<class Type>
void Foam::lookupTable3D<Type>::setY
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
void Foam::lookupTable3D<Type>::setY
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
}


template<class Type>
void Foam::lookupTable3D<Type>::setZ
(
    const Field<scalar>& z,
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
    const Field<scalar>& z,
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
}


template<class Type>
void Foam::lookupTable3D<Type>::setData
(
    const Field<Field<Field<Type>>>& data,
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

    realDataPtr_ = new Field<Field<Field<Type>>>(data);
    data_ = data;

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
void Foam::lookupTable3D<Type>::setData
(
    const Field<Field<Field<Type>>>& data,
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
    ijk_.x() = xIndexing_->findIndex(modX_()(x), xModValues_);
    ijk_.y() = yIndexing_->findIndex(modY_()(y), yModValues_);
    ijk_.z() = zIndexing_->findIndex(modZ_()(z), zModValues_);
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

    ijk_.x() = xIndexing_->findIndex(xMod, xModValues_);
    ijk_.y() = yIndexing_->findIndex(yMod, yModValues_);
    ijk_.z() = zIndexing_->findIndex(zMod, zModValues_);

    labelList is, js, ks;
    scalarList wxs, wys, wzs;

    xInterpolator_->updateWeights(xMod, ijk_.x(), xModValues_, is, wxs);
    yInterpolator_->updateWeights(yMod, ijk_.y(), yModValues_, js, wys);
    zInterpolator_->updateWeights(zMod, ijk_.z(), zModValues_, ks, wzs);

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
        xModValues_.size()
    );

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
        yModValues_.size()
    );

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
        zModValues_.size()
    );

    const dictionary& fDict(dict.optionalSubDict(name + "Coeffs"));

    word modType = fDict.lookupOrDefault<word>("mod", "none");
    isReal = fDict.lookupOrDefault<Switch>("isReal", true);

    fileName file(fDict.lookup<word>("file"));
    Field<Field<Field<Type>>> data
    (
        xModValues_.size(),
        Field<Field<Type>>
        (
            yModValues_.size(),
            Field<Type>(zModValues_.size())
        )
    );
    read3DTable
    (
        file,
        dict.lookupOrDefault<string>("delim", ","),
        dict.lookupOrDefault<string>("rowDelim", ";"),
        data,
        dict.lookupOrDefault<Switch>("flipTable", true),
        !canRead
    );
    setData(data, modType, isReal);
}

// ************************************************************************* //
