/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
#include "DynamicList.H"
#include "Field.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable2D<Type>::readTable
(
    const fileName& file,
    const string& delim,
    Field<Field<Type>>& data
)
{
    fileName fNameExpanded(file);
    fNameExpanded.expand();

    // Open a stream and check it
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fNameExpanded));
    ISstream& is = isPtr();
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file" << file << nl
            << exit(FatalIOError);
    }

    DynamicList<Tuple2<scalar, scalar>> values;

    label i = 0;
    while (is.good())
    {
        string line;
        is.getLine(line);

        string lineEntry = line;
        lineEntry.replaceAll(delim, " ");
        lineEntry = '(' + lineEntry + ')';
        IStringStream iss(lineEntry);

        List<Type> vals(iss);

        if (vals.size() <= 1)
        {
            break;
        }

        if (vals.size() != yValues_.size())
        {
            FatalErrorInFunction
                << "Incompatible table rows" << endl
                << line
                << abort(FatalError);
        }
        forAll(vals, j)
        {
            data[i][j] = vals[j];
        }
        i++;
    }
}


template<class Type>
void Foam::lookupTable2D<Type>::findUniformIndexes
(
    const scalar xy,
    const scalarField& XY,
    label& IJ,
    scalar& f
)
{
    scalar xyMin = XY[0];
    scalar dxy = XY[1] - XY[0];

    scalar nxy = XY.size();
    scalar ij = (xy - xyMin)/dxy;
    if (ij < 0)
    {
        IJ = 0;
        f = 1.0;
        return;
    }
    IJ = floor(ij);
    if (IJ >= nxy - 1)
    {
        IJ = nxy - 2;
        f = 0.0;
        return;
    }

    f = 1.0 - (ij - floor(IJ));
    return;
}


template<class Type>
void Foam::lookupTable2D<Type>::findNonuniformIndexes
(
    const scalar xy,
    const scalarField& XY,
    label& IJ,
    scalar& f
)
{
    if (xy < XY[0])
    {
        IJ = 0;
        f = 1.0;
        return;
    }

    forAll(XY, ij)
    {
        if (xy < XY[ij] && xy < XY[ij+1])
        {
            IJ = ij;
            f = 1.0 - (xy - XY[ij])/(XY[ij+1] - XY[ij]);
            return;
        }
    }
    IJ = XY.size() - 2;
    f = 0.0;
    return;
}


template<class Type>
Type Foam::lookupTable2D<Type>::getValue
(
    const label ij,
    const scalar f,
    const Field<Type>& xy
) const
{
    if (ij >= xy.size())
    {
        return xy.last();
    }

    return xy[ij] + f*(xy[ij+1] - xy[ij]);
}


template<class Type>
bool Foam::lookupTable2D<Type>::checkUniform(const scalarField& xy) const
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
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    xValues_(),
    yValues_(),
    findXIndex_(nullptr),
    findYIndex_(nullptr)
{}

template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    data_(),
    xModValues_(),
    yModValues_(),
    xValues_(),
    yValues_(),
    findXIndex_(nullptr),
    findYIndex_(nullptr)
{
    read(dict, xName, yName, name);
}


template<class Type>
Foam::lookupTable2D<Type>::lookupTable2D
(
    const scalarField& x,
    const scalarField& y,
    const Field<Field<Type>>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const bool isReal
)
:
    modFunc_(nullptr),
    invModFunc_(nullptr),
    modXFunc_(nullptr),
    invModXFunc_(nullptr),
    modYFunc_(nullptr),
    invModYFunc_(nullptr),
    data_(data),
    xModValues_(x),
    yModValues_(y),
    xValues_(x),
    yValues_(y),
    findXIndex_(nullptr),
    findYIndex_(nullptr)
{
    set(x, y, data, modXType, modYType, modType, isReal);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::lookupTable2D<Type>::~lookupTable2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lookupTable2D<Type>::set
(
    const scalarField& x,
    const scalarField& y,
    const Field<Field<Type>>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const bool isReal
)
{
    setX(x, modXType, isReal);
    setY(y, modYType, isReal);
    setData(data, modType, isReal);
}


template<class Type>
void Foam::lookupTable2D<Type>::setX
(
    const scalarField& x,
    const word& modXType,
    const bool isReal
)
{
    setMod(modXType, modXFunc_, invModXFunc_);

    if (isReal)
    {
        forAll(xValues_, i)
        {
            xValues_[i] = invModXFunc_(xModValues_[i]);
        }
    }
    else
    {
        forAll(xValues_, i)
        {
            xModValues_[i] = modXFunc_(xValues_[i]);
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
    const scalarField& y,
    const word& modYType,
    const bool isReal
)
{
    setMod(modYType, modYFunc_, invModYFunc_);

    if (isReal)
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
    const word& modType,
    const bool isReal
)
{
    setMod(modType, modFunc_, invModFunc_);
    data_ = data;

    if (!isReal)
    {
        forAll(data_, i)
        {
            forAll(data_[i], j)
            {
                data_[i][j] = modFunc_(data_[i][j]);
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Foam::Field<Type>>>
Foam::lookupTable2D<Type>::realData() const
{
    tmp<Field<Field<Type>>> tmpf(new Field<Field<Type>>(data_));
    Field<Field<Type>>& f = tmpf.ref();
    forAll(f, i)
    {
        forAll(f[i], j)
        {
            f[i][j] = invModFunc_(f[i][j]);
        }
    }
    return tmpf;
}


template<class Type>
void Foam::lookupTable2D<Type>::update(const scalar x, const scalar y) const
{
    findXIndex_(modXFunc_(x), xModValues_, i_, fx_);
    findYIndex_(modYFunc_(y), yModValues_, j_, fy_);
}


template<class Type>
void Foam::lookupTable2D<Type>::updateX(const scalar x) const
{
    findXIndex_(modXFunc_(x), xModValues_, i_, fx_);
}


template<class Type>
void Foam::lookupTable2D<Type>::updateY(const scalar y) const
{
    findYIndex_(modYFunc_(y), yModValues_, j_, fy_);
}


template<class Type>
Type Foam::lookupTable2D<Type>::lookup(const scalar x, const scalar y) const
{
    update(x, y);
    return
        invModFunc_
        (
            data_[i_][j_]*fx_*fy_
          + data_[i_+1][j_]*(1.0 - fx_)*fy_
          + data_[i_][j_+1]*fx_*(1.0 - fy_)
          + data_[i_+1][j_+1]*(1.0 - fx_)*(1.0 - fy_)
        );
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdX(const scalar x, const scalar y) const
{
    update(x, y);

    const Type& mm(data_[i_][j_]);
    const Type& pm(data_[i_+1][j_]);
    const Type& mp(data_[i_][j_+1]);
    const Type& pp(data_[i_+1][j_+1]);

    return
        (
            invModFunc_(pm*fy_ + pp*(1.0 - fy_))
          - invModFunc_(mm*fy_ + mp*(1.0 - fy_))
        )/(xValues_[i_+1] - xValues_[i_]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::dFdY(const scalar x, const scalar y) const
{
    update(x, y);

    const Type& mm(data_[i_][j_]);
    const Type& pm(data_[i_+1][j_]);
    const Type& mp(data_[i_][j_+1]);
    const Type& pp(data_[i_+1][j_+1]);

    return
        (
            invModFunc_(mp*fx_ + pp*(1.0 - fx_))
          - invModFunc_(mm*fx_ + pm*(1.0 - fx_))
        )/(yValues_[j_+1] - yValues_[j_]);
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdX2(const scalar x, const scalar y) const
{
    update(x, y);

    if (i_ == 0)
    {
        i_++;
    }

    const Type gmm(invModFunc_(data_[i_-1][j_]));
    const Type gm(invModFunc_(data_[i_][j_]));
    const Type gpm(invModFunc_(data_[i_+1][j_]));

    const Type gmp(invModFunc_(data_[i_-1][j_+1]));
    const Type gp(invModFunc_(data_[i_][j_+1]));
    const Type gpp(invModFunc_(data_[i_+1][j_+1]));

    const scalar& xm(xValues_[i_-1]);
    const scalar& xi(xValues_[i_]);
    const scalar& xp(xValues_[i_+1]);

    const Type gPrimepm((gpm - gm)/(xp - xi));
    const Type gPrimemm((gm - gmm)/(xi - xm));
    const Type gPrimepp((gpp - gp)/(xp - xi));
    const Type gPrimemp((gp - gmp)/(xi - xm));

    return
        fy_*(gPrimepm - gPrimemm)/(0.5*(xp - xm))
      + (1.0 - fy_)*(gPrimepp - gPrimemp)/(0.5*(xp - xm));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdY2(const scalar x, const scalar y) const
{
    update(x, y);

    if (j_ == 0)
    {
        j_++;
    }

    const Type gmm(invModFunc_(data_[i_][j_-1]));
    const Type gm(invModFunc_(data_[i_][j_]));
    const Type gmp(invModFunc_(data_[i_][j_+1]));

    const Type gpm(invModFunc_(data_[i_+1][j_-1]));
    const Type gp(invModFunc_(data_[i_+1][j_]));
    const Type gpp(invModFunc_(data_[i_+1][j_+1]));

    const scalar& ym(yValues_[j_-1]);
    const scalar& yi(yValues_[j_]);
    const scalar& yp(yValues_[j_+1]);

    const Type gPrimemp((gmp - gm)/(yp - yi));
    const Type gPrimemm((gm - gmm)/(yi - ym));
    const Type gPrimepp((gpp - gp)/(yp - yi));
    const Type gPrimepm((gp - gpm)/(yi - ym));

    return
        fx_*(gPrimemp - gPrimemm)/(0.5*(yp - ym))
      + (1.0 - fx_)*(gPrimepp - gPrimepm)/(0.5*(yp - ym));
}


template<class Type>
Type Foam::lookupTable2D<Type>::d2FdXdY(const scalar x, const scalar y) const
{
    update(x, y);

    const Type gmm(invModFunc_(data_[i_][j_]));
    const Type gmp(invModFunc_(data_[i_][j_+1]));
    const Type gpm(invModFunc_(data_[i_+1][j_]));
    const Type gpp(invModFunc_(data_[i_+1][j_+1]));

    const scalar& xm(xValues_[i]);
    const scalar& xp(xValues_[i+1]);

    const scalar& ym(yValues_[j]);
    const scalar& yp(yValues_[j+1]);

    return ((gpp - gmp)/(xp - xm) - (gpm - gmm)/(xp - xm))/(yp - ym);
}


template<class Type>
void Foam::lookupTable2D<Type>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name
)
{
    setMod(dict.lookup("mod"), modFunc_, invModFunc_);

    Switch isReal(dict.lookupOrDefault<Switch>("isReal", true));

    if (dict.found(xName))
    {
        xValues_ = dict.lookup<scalarList>(xName);
        setMod(dict.lookup(xName + "Mod"), modXFunc_, invModXFunc_);
    }
    else
    {
        const dictionary& xDict(dict.subDict(xName + "Coeffs"));
        label nx = xDict.lookup<label>("n");
        scalar dx = xDict.lookup<scalar>("delta");
        scalar minx = xDict.lookup<scalar>("min");

        setMod(xDict.lookup("mod"), modXFunc_, invModXFunc_);

        xValues_.resize(nx);
        forAll(xValues_, i)
        {
            xValues_[i] = minx + dx*i;
        }
    }

    {
        if (dict.found(yName))
        {
            yValues_ = dict.lookup<scalarList>(yName);
            yModValues_.resize(yValues_.size());

            setMod(dict.lookup(yName + "Mod"), modYFunc_, invModYFunc_);
        }
        else
        {
            const dictionary& yDict(dict.subDict(yName + "Coeffs"));
            label ny = yDict.lookup<label>("n");
            scalar dy = yDict.lookup<scalar>("delta");
            scalar miny = yDict.lookup<scalar>("min");

            setMod(dict.lookup("mod"), modYFunc_, invModYFunc_);

            yValues_.resize(ny);
            forAll(yValues_, j)
            {
                yValues_[j] = miny + dy*j;
            }
        }
    }

    data_.resize(xValues_.size());
    forAll(data_, i)
    {
        data_[i].resize(yValues_.size());
    }


    fileName file(dict.lookup<word>("file"));
    readTable
    (
        file,
        dict.lookupOrDefault<string>("delim", ","),
        data_
    );

    if (!isReal)
    {
        xModValues_ = xValues_;
        yModValues_ = yValues_;
        forAll(xValues_, i)
        {
            xValues_[i] = invModYFunc_(xValues_[i]);
        }
        forAll(yValues_, i)
        {
            yValues_[i] = invModYFunc_(yValues_[i]);
        }
    }
    else
    {
        xModValues_.resize(xValues_.size());
        yModValues_.resize(yValues_.size());
        forAll(xValues_, i)
        {
            xModValues_[i] = modYFunc_(xValues_[i]);
            forAll(yValues_, j)
            {
                data_[i][j] = modFunc_(data_[i][j]);
            }
        }
        forAll(yValues_, i)
        {
            yModValues_[i] = modYFunc_(yValues_[i]);
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

    if (checkUniform(yModValues_))
    {
        findYIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        findYIndex_ = &lookupTable2D::findNonuniformIndexes;
    }
}

// ************************************************************************* //
