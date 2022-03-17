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


template<class Type>
bool Foam::lookupTable3D<Type>::checkUniform(const List<scalar>& xyz) const
{
    scalar dxyz = xyz[1] - xyz[0];

    for (label i = 2; i < xyz.size(); i++)
    {
        if (mag(xyz[i] - xyz[i-1] - dxyz) > small)
        {
            return false;
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable3D<Type>::lookupTable3D()
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
    modZType_("none"),
    needZMod_(false),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
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
    modZType_("none"),
    needZMod_(false),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
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
    modZType_("none"),
    needZMod_(false),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpType_("linearClamp"),
    interpFunc_(nullptr),
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
        interpolationScheme,
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
    const word& interpolationScheme,
    const bool isReal
)
{
    setX(x, modXType, isReal);
    setY(y, modYType, isReal);
    setZ(z, modZType, isReal);
    setData(data, modType, isReal);

    interpType_ = interpolationScheme;
    setInterp(interpolationScheme, interpFunc_);
}


template<class Type>
void Foam::lookupTable3D<Type>::setX
(
    const Field<scalar>& x,
    const word& modX,
    const bool isReal
)
{
    modXType_ = modX;
    needXMod_ = modX != "none";

    setMod<scalar>(modX, modXFunc_, invModXFunc_);
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
        findXIndex_ = &lookupTable3D::findUniformIndexes;
    }
    else
    {
        findXIndex_ = &lookupTable3D::findNonuniformIndexes;
    }
}


template<class Type>
void Foam::lookupTable3D<Type>::setY
(
    const Field<scalar>& y,
    const word& modY,
    const bool isReal
)
{
    modYType_ = modY;
    needYMod_ = modY != "none";

    setMod<scalar>(modY, modYFunc_, invModYFunc_);
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
        findYIndex_ = &lookupTable3D::findUniformIndexes;
    }
    else
    {
        findYIndex_ = &lookupTable3D::findNonuniformIndexes;
    }
}


template<class Type>
void Foam::lookupTable3D<Type>::setZ
(
    const Field<scalar>& z,
    const word& modZ,
    const bool isReal
)
{
    modZType_ = modZ;
    needZMod_ = modZ != "none";

    setMod<scalar>(modZ, modZFunc_, invModZFunc_);
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

    if (!needZMod_)
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
                zModValues_[i] = modZFunc_(z[i]);
            }
        }
        else
        {
            forAll(z, i)
            {
                (*zValuesPtr_)[i] = invModZFunc_(z[i]);
            }
        }
    }

    if (checkUniform(zModValues_))
    {
        findZIndex_ = &lookupTable3D::findUniformIndexes;
    }
    else
    {
        findZIndex_ = &lookupTable3D::findNonuniformIndexes;
    }
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

    if (!needMod_)
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
                    data_[i][j][k] = modFunc_(data[i][j][k]);
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
                    (*realDataPtr_)[i][j][k] = invModFunc_(data[i][j][k]);
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
    modType_ = mod;
    needMod_ = mod != "none";

    setMod<Type>(mod, modFunc_, invModFunc_);
    setData(data, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::update
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    scalar xMod(modXFunc_(x));
    scalar yMod(modYFunc_(y));
    scalar zMod(modZFunc_(z));

    ijk_.x() = findXIndex_(xMod, xModValues_);
    ijk_.y() = findYIndex_(yMod, yModValues_);
    ijk_.z() = findZIndex_(zMod, zModValues_);

    labelList is, js, ks;
    scalarList wxs, wys, wzs;

    interpFunc_(x, ijk_.x(), xModValues_, is, wxs);
    interpFunc_(y, ijk_.y(), yModValues_, js, wys);
    interpFunc_(z, ijk_.z(), zModValues_, ks, wzs);

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
    return invModFunc_(modf);
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
    interpType_ =
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp");
    setInterp(interpType_, interpFunc_);

    scalarField x, y, z;
    bool isReal = readComponent
    (
        dict,
        xName,
        modXType_,
        x
    );
    setX(x, modXType_, isReal);

    isReal = readComponent
    (
        dict,
        yName,
        modYType_,
        y
    );
    setY(y, modYType_, isReal);

    isReal = readComponent
    (
        dict,
        zName,
        modZType_,
        z
    );
    setZ(z, modZType_, isReal);

    const dictionary& fDict(dict.optionalSubDict(name + "Coeffs"));

    modType_ = fDict.lookupOrDefault<word>("mod", "none");
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
    setData(data, modZType_, isReal);
}

// ************************************************************************* //
