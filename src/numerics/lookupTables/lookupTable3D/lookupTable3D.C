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
    modType_("none"),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_("none"),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYType_("none"),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    modZType_("none"),
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
    const word& name,
    const bool canRead
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
    interpType_("linearClamp"),
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
    modType_(modType),
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXType_(modXType),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYType_(modYType),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    modZType_(modZType),
    modZFunc_(nullptr),
    invModZFunc_(nullptr),
    findXIndex_(nullptr),
    findYIndex_(nullptr),
    findZIndex_(nullptr),
    interpType_(interpolationScheme),
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

    interpType_ = interpolationScheme;
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
    modXType_ = modXType;
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
    const Field<scalar>& y,
    const word& modYType,
    const bool isReal
)
{
    modYType_ = modYType;
    setMod(modYType, modYFunc_, invModYFunc_);
    setY(y, isReal);
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
    modZType_ = modZType;
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
    modType_ = modType;
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
    const word& name,
    const bool canRead
)
{
    interpType_ =
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp");
    setInterp(interpType_, interpFunc_);

    const dictionary& fDict(dict.subDict(name + "Coeffs"));

    modType_ = fDict.lookupOrDefault<word>("mod", "none");
    setMod(modType_, modFunc_, invModFunc_);

    Switch isReal(fDict.lookupOrDefault<Switch>("isReal", true));

    readComponent
    (
        dict,
        xName,
        modXType_,
        xValues_,
        xModValues_,
        modXFunc_,
        invModXFunc_,
        findXIndex_,
        canRead
    );
    readComponent
    (
        dict,
        yName,
        modYType_,
        yValues_,
        yModValues_,
        modYFunc_,
        invModYFunc_,
        findYIndex_,
        canRead
    );
    readComponent
    (
        dict,
        zName,
        modZType_,
        zValues_,
        zModValues_,
        modZFunc_,
        invModZFunc_,
        findZIndex_,
        canRead
    );

    if (canRead)
    {
        data_.resize(xValues_.size());
        forAll(data_, i)
        {
            data_[i].resize(yValues_.size());
            forAll(data_[i], j)
            {
                data_[i][j].resize(zValues_.size());
            }
        }
    }


    fileName file(fDict.lookup<word>("file"));
    read3DTable
    (
        file,
        dict.lookupOrDefault<string>("delim", ","),
        dict.lookupOrDefault<string>("rowDelim", ";"),
        data_,
        dict.lookupOrDefault<Switch>("flipTable", true),
        !canRead
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
    word& type,
    Field<scalar>& values,
    Field<scalar>& modValues,
    modFuncType& modFunc,
    modFuncType& invModFunc,
    findIndexFunc& findIndex,
    const bool canRead
)
{
    Switch isReal = true;
    bool canSetMod = true;
    if (parentDict.found(name + "Coeffs"))
    {
        const dictionary& dict(parentDict.subDict(name + "Coeffs"));
        type = dict.lookupOrDefault<word>("mod", "none");
        setMod
        (
            type,
            modFunc,
            invModFunc
        );

        if (dict.found(name))
        {
            if (canRead)
            {
                values = dict.lookup<Field<Type>>(name);
                modValues.resize(values.size());
                isReal = dict.lookupOrDefault<Switch>("isReal", true);
            }
            else
            {
                canSetMod = false;
            }
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
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
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
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
        }
    }
    else if (parentDict.found(name))
    {
        if (canRead)
        {
            values = parentDict.lookup<Field<Type>>(name);
            modValues.resize(values.size());
            isReal =
                parentDict.lookupOrDefault<Switch>
                (
                    name + "isReal",
                    true
                );
        }
        else
        {
            canSetMod = false;
        }
        type = parentDict.lookupOrDefault<word>(name + "Mod", "none");
        setMod
        (
            type,
            modFunc,
            invModFunc
        );
    }
    else
    {
        FatalErrorInFunction
            << "Either a list of values or " << name << "Coeffs must" << nl
            << "be provided for " << name << endl
            << abort(FatalError);
    }

    if (canSetMod)
    {
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
}

// ************************************************************************* //
