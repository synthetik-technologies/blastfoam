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
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xValues_(),
    yValues_(),
    zValues_(),
    i_(0),
    j_(0),
    k_(0),
    fx_(0),
    fy_(0),
    fz_(0)
{}


template<class Type>
Foam::lookupTable3D<Type>::lookupTable3D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& zName,
    const word& name
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xValues_(),
    yValues_(),
    zValues_(),
    i_(0),
    j_(0),
    k_(0),
    fx_(0),
    fy_(0),
    fz_(0)
{
    read(dict, xName, yName, zName, name);
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
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    zModValues_(),
    xValues_(),
    yValues_(),
    zValues_(),
    i_(0),
    j_(0),
    k_(0),
    fx_(0),
    fy_(0),
    fz_(0)
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
{}


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
    setInterp(interpolationScheme, interpFunc_);
}


template<class Type>
void Foam::lookupTable3D<Type>::setX
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
void Foam::lookupTable3D<Type>::setX
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
    const Field<scalar>& x,
    const word& modXType,
    const bool isReal
)
{
    setMod(modXType, modXFunc_, invModXFunc_);
    setX(x, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::setY
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
    const word& modZType,
    const bool isReal
)
{
    setMod(modZType, modZFunc_, invModZFunc_);
    setZ(z, isReal);
}


template<class Type>
void Foam::lookupTable3D<Type>::setZ
(
    const Field<scalar>& z,
    const bool isReal
)
{
    if (!isReal)
    {
        forAll(zValues_, j)
        {
            zValues_[j] = invModZFunc_(zModValues_[j]);
        }
    }
    else
    {
        forAll(zValues_, j)
        {
            zModValues_[j] = modZFunc_(zValues_[j]);
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
    data_.resize(data.size());
    forAll(data, i)
    {
        data_[i].resize(data[i].size());
        forAll(data[i], j)
        {
            data_[i][j] = data[i][j];
        }
    }

    if (isReal)
    {
        forAll(data_, i)
        {
            forAll(data_[i], j)
            {
                forAll(data_[i][j], k)
                {
                    data_[i][j][k] = modFunc_(data_[i][j][k]);
                }
            }
        }
    }
}


template<class Type>
void Foam::lookupTable3D<Type>::setData
(
    const Field<Field<Field<Type>>>& data,
    const word& modType,
    const bool isReal
)
{
    setMod(modType, modFunc_, invModFunc_);
    setData(data, isReal);
}


template<class Type>
Foam::tmp<Foam::Field<Foam::Field<Foam::Field<Type>>>>
Foam::lookupTable3D<Type>::realData() const
{
    tmp<Field<Field<Field<Type>>>> tmpf
    (
        new Field<Field<Field<Type>>>(data_)
    );
    Field<Field<Field<Type>>>& f = tmpf.ref();
    forAll(f, i)
    {
        forAll(f[i], j)
        {
            forAll(f[i][j], k)
            {
                f[i][j][k] = invModFunc_(f[i][j][k]);
            }
        }
    }
    return tmpf;
}


template<class Type>
void Foam::lookupTable3D<Type>::update
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    updateX(x);
    updateY(y);
    updateZ(z);
}


template<class Type>
void Foam::lookupTable3D<Type>::updateX(const scalar x) const
{
    i_ = findXIndex_(modXFunc_(x), xModValues_);
}


template<class Type>
void Foam::lookupTable3D<Type>::updateY(const scalar y) const
{
    j_ = findYIndex_(modYFunc_(y), yModValues_);

}


template<class Type>
void Foam::lookupTable3D<Type>::updateZ(const scalar z) const
{
    k_ = findZIndex_(modZFunc_(z), zModValues_);

}


template<class Type>
void Foam::lookupTable3D<Type>::updateWeights
(
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    updateXWeight(x);
    updateYWeight(y);
    updateZWeight(z);
}


template<class Type>
void Foam::lookupTable3D<Type>::updateXWeight(const scalar x) const
{
    fx_ = linearWeight(modXFunc_(x), xModValues_[i_], xModValues_[i_ + 1]);
}


template<class Type>
void Foam::lookupTable3D<Type>::updateYWeight(const scalar y) const
{
    fy_ = linearWeight(modYFunc_(y), yModValues_[j_], yModValues_[j_ + 1]);

}


template<class Type>
void Foam::lookupTable3D<Type>::updateZWeight(const scalar z) const
{
    fz_ = linearWeight(modZFunc_(z), zModValues_[k_], zModValues_[k_ + 1]);

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
    return
        invModFunc_
        (
            interpFunc_
            (
                modXFunc_(x), modYFunc_(y), modZFunc_(z),
                i_, j_, k_,
                xModValues_, yModValues_, zModValues_,
                data_
            )
        );
}


template<class Type>
void Foam::lookupTable3D<Type>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& zName,
    const word& name
)
{
    word interpolationScheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );
    setInterp(interpolationScheme, interpFunc_);

    const dictionary& fDict(dict.subDict(name + "Coeffs"));
    setMod(fDict.lookup("mod"), modFunc_, invModFunc_);
    Switch isReal(fDict.lookupOrDefault<Switch>("isReal", true));

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
    readComponent
    (
        dict,
        zName,
        zValues_,
        zModValues_,
        modZFunc_,
        invModZFunc_,
        findZIndex_
    );

    data_.resize(xValues_.size());
    forAll(data_, i)
    {
        data_[i].resize(yValues_.size());
        forAll(data_[i], j)
        {
            data_[i][j].resize(zValues_.size());
        }
    }


    fileName file(fDict.lookup<word>("file"));
    read3DTable
    (
        file,
        dict.lookupOrDefault<string>("delim", ","),
        dict.lookupOrDefault<string>("rowDelim", ";"),
        data_,
        dict.lookupOrDefault<Switch>("flipTable", true)
    );
    if (isReal)
    {
        forAll(xValues_, i)
        {
            forAll(yValues_, j)
            {
                forAll(zValues_, k)
                {
                    data_[i][j][k] = modFunc_(data_[i][j][k]);
                }
            }
        }
    }
}


template<class Type>
void Foam::lookupTable3D<Type>::readComponent
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
    const dictionary& dict(parentDict.subDict(name + "Coeffs"));
    Switch isReal(dict.lookupOrDefault<Switch>("isReal", true));
    setMod(dict.lookupOrDefault<word>("mod", "none"), modFunc, invModFunc);

    if (dict.found(name))
    {
        values = dict.lookup<Field<scalar>>(name);
    }
    else if (dict.found("file"))
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
        scalar dy = dict.lookup<scalar>("delta");
        scalar miny = dict.lookup<scalar>("min");

        values.resize(ny);
        forAll(values, j)
        {
            values[j] = miny + dy*j;
        }
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
        findIndex = &lookupTable3D::findUniformIndexes;
    }
    else
    {
        findIndex = &lookupTable3D::findNonuniformIndexes;
    }
}

// ************************************************************************* //
