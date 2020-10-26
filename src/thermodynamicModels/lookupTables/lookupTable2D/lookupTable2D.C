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

Foam::scalar Foam::lookupTable2D::readValue(const List<string>& split) const
{
    if (0 >= split.size())
    {
        FatalErrorInFunction
            << "No column " << 0 << " in "
            << split << endl
            << exit(FatalError);
    }

    return readScalar(IStringStream(split[0])());
}


void Foam::lookupTable2D::readTable(const fileName& file, Field<scalarField>& data)
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

    bool mergeSeparators_=false;
    string separator_=';';

    label i = 0;
    while (is.good())
    {
        string line;
        is.getLine(line);

        label n = 0;
        std::size_t pos = 0;
        DynamicList<string> split;

        if (mergeSeparators_)
        {
            std::size_t nPos = 0;

            while ((pos != std::string::npos) && (n <= ny_))
            {
                bool found = false;
                while (!found)
                {
                    nPos = line.find(separator_, pos);

                    if ((nPos != std::string::npos) && (nPos - pos == 0))
                    {
                        pos = nPos + 1;
                    }
                    else
                    {
                        found = true;
                    }
                }

                nPos = line.find(separator_, pos);

                if (nPos == std::string::npos)
                {
                    split.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    split.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }
        else
        {
            while ((pos != std::string::npos) && (n <= ny_))
            {
                std::size_t nPos = line.find(separator_, pos);

                if (nPos == std::string::npos)
                {
                    split.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    split.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }

        if (split.size() <= 1)
        {
            break;
        }
        for (label j = 0; j < ny_; j++)
        {
            data[i][j] = readScalar(IStringStream(split[j])());
        }
        i++;
    }
}

void Foam::lookupTable2D::findUniformIndexes
(
    const scalar& xy,
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


void Foam::lookupTable2D::findNonuniformIndexes
(
    const scalar& xy,
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


Foam::label Foam::lookupTable2D::boundi
(
    const scalar& f,
    const Field<scalarField>& data,
    const label j
) const
{
    for (label i = 0; i < data.size() - 1; i++)
    {
        if (f > data[i][j] && f < data[i+1][j])
        {
            return i;
        }
    }
    return data.size() - 2;
}


Foam::label Foam::lookupTable2D::boundj
(
    const scalar& f,
    const Field<scalarField>& data,
    const label i
) const
{
    for (label j = 0; j < data[i].size() - 1; j++)
    {
        if (f > data[i][j] && f < data[i][j+1])
        {
            return j;
        }
    }
    return data[i].size() - 2;
}


void Foam::lookupTable2D::boundij
(
    const scalar& f,
    const Field<scalarField>& data,
    label& i,
    label& j
) const
{
    for (i = 0; i < data.size() - 1; i++)
    {
        for (j = 0; j < data[i].size() - 1; j++)
        {
            if
            (
                f > data[i][j]
             && f < data[i+1][j]
             && f < data[i][j+1]
             && f < data[i+1][j+1]
            )
            {
                return;
            }
        }
    }
    i--;
    j--;
    return;
}


Foam::scalar Foam::lookupTable2D::getValue
(
    const label ij,
    const scalar& f,
    const scalarField& xy
) const
{
    if (ij >= xy.size())
    {
        return xy.last();
    }

    return xy[ij] + f*(xy[ij+1] - xy[ij]);
}


bool Foam::lookupTable2D::checkUniform(const scalarField& xy) const
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

Foam::lookupTable2D::lookupTable2D
(
    const dictionary& dict,
    const word& xName,
    const word& yName
)
:
    modType_(dict.lookup("mod")),
    modFunc_(NULL),
    invModFunc_(NULL),
    modXType_(dict.lookup(xName + "Mod")),
    modXFunc_(NULL),
    invModXFunc_(NULL),
    modYType_(dict.lookup(yName + "Mod")),
    modYFunc_(NULL),
    invModYFunc_(NULL),
    uniformX_(dict.lookupOrDefault(xName + "Uniform", true)),
    uniformY_(dict.lookupOrDefault(yName + "Uniform", true)),
    nx_(dict.lookupType<label>("n" + xName)),
    ny_(dict.lookupType<label>("n" + yName)),
    data_(nx_, scalarField(ny_, 0.0)),
    xMod_(nx_, 0.0),
    yMod_(ny_, 0.0),
    x_(nx_, 0.0),
    y_(ny_, 0.0),
    findXIndex_(NULL),
    findYIndex_(NULL)
{
    setMod(modType_, modFunc_, invModFunc_);
    setMod(modXType_, modXFunc_, invModXFunc_);
    setMod(modYType_, modYFunc_, invModYFunc_);

    fileName file(dict.lookupType<word>("file"));
    readTable(file, data_);

    if (uniformX_)
    {
        scalar dx(dict.lookupType<scalar>("d" + xName));
        scalar minx(dict.lookupType<scalar>("min" + xName));
        forAll(xMod_, i)
        {
            xMod_[i] = minx + i*dx;
        }
        findXIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        xMod_ = dict.lookupType<scalarField>(xName);
        findXIndex_ = &lookupTable2D::findNonuniformIndexes;
    }

    if (uniformY_)
    {
        scalar dy(dict.lookupType<scalar>("d" + yName));
        scalar miny(dict.lookupType<scalar>("min" + yName));
        forAll(yMod_, j)
        {
            yMod_[j] = miny + dy*j;
        }
        findYIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        yMod_ = dict.lookupType<scalarField>(yName);
        findYIndex_ = &lookupTable2D::findNonuniformIndexes;
    }

    forAll(x_, i)
    {
        x_[i] = invModXFunc_(xMod_[i]);
    }
    forAll(y_, j)
    {
        y_[j] = invModYFunc_(yMod_[j]);
    }
}

Foam::lookupTable2D::lookupTable2D
(
    const fileName& file,
    const word& modType,
    const word& modXType,
    const word& modYType,
    const label nx,
    const label ny,
    const scalar& xMin,
    const scalar& dx,
    const scalar& yMin,
    const scalar& dy
)
:
    modType_(modType),
    modFunc_(NULL),
    invModFunc_(NULL),
    modXType_(modXType),
    modXFunc_(NULL),
    invModXFunc_(NULL),
    modYType_(modYType),
    modYFunc_(NULL),
    invModYFunc_(NULL),
    uniformX_(true),
    uniformY_(true),
    nx_(nx),
    ny_(ny),
    data_(nx_, scalarField(ny_, 0.0)),
    xMod_(nx_, 0.0),
    yMod_(ny_, 0.0),
    x_(nx_, 0.0),
    y_(ny_, 0.0),
    findXIndex_(NULL),
    findYIndex_(NULL)
{
    setMod(modType, modFunc_, invModFunc_);
    setMod(modXType, modXFunc_, invModXFunc_);
    setMod(modYType, modYFunc_, invModYFunc_);

    readTable(file, data_);

    forAll(xMod_, i)
    {
        xMod_[i] = xMin + i*dx;
        x_[i] = invModXFunc_(xMod_[i]);
    }
    findXIndex_ = &lookupTable2D::findUniformIndexes;

    forAll(yMod_, j)
    {
        yMod_[j] = yMin + dy*j;
        y_[j] = invModYFunc_(yMod_[j]);
    }
    findYIndex_ = &lookupTable2D::findUniformIndexes;
}


Foam::lookupTable2D::lookupTable2D
(
    const Field<scalarField>& data,
    const scalarField& x,
    const scalarField& y,
    const word& modType,
    const word& modXType,
    const word& modYType,
    const bool modified
)
:
    modType_(modType),
    modFunc_(NULL),
    invModFunc_(NULL),
    modXType_(modXType),
    modXFunc_(NULL),
    invModXFunc_(NULL),
    modYType_(modYType),
    modYFunc_(NULL),
    invModYFunc_(NULL),
    nx_(data.size()),
    ny_(data[0].size()),
    data_(data),
    xMod_(modified ? x : scalarField(nx_, 0)),
    yMod_(modified ? y : scalarField(ny_, 0)),
    x_(!modified ? x : scalarField(nx_, 0)),
    y_(!modified ? y : scalarField(ny_, 0)),
    findXIndex_(NULL),
    findYIndex_(NULL)
{
    setMod(modType, modFunc_, invModFunc_);
    setMod(modXType, modXFunc_, invModXFunc_);
    setMod(modYType, modYFunc_, invModYFunc_);

    if (modified)
    {
        forAll(x_, i)
        {
            x_[i] = invModXFunc_(xMod_[i]);
        }
        forAll(y_, j)
        {
            y_[j] = invModYFunc_(yMod_[j]);
        }
    }
    else
    {
        forAll(x_, i)
        {
            xMod_[i] = modXFunc_(x_[i]);
        }
        forAll(y_, j)
        {
            yMod_[j] = modYFunc_(y_[j]);
        }

        forAll(data_, i)
        {
            forAll(data_[i], j)
            {
                data_[i][j] = modFunc_(data_[i][j]);
            }
        }
    }

    uniformX_ = checkUniform(xMod_);
    uniformY_ = checkUniform(yMod_);

    if (uniformX_)
    {
        findXIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        findXIndex_ = &lookupTable2D::findNonuniformIndexes;
    }

    if (uniformY_)
    {
        findYIndex_ = &lookupTable2D::findUniformIndexes;
    }
    else
    {
        findYIndex_ = &lookupTable2D::findNonuniformIndexes;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lookupTable2D::~lookupTable2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalarField>> Foam::lookupTable2D::realData() const
{
    tmp<Field<scalarField>> tmpf(new Field<scalarField>(data_));
    Field<scalarField>& f = tmpf.ref();
    forAll(f, i)
    {
        forAll(f[i], j)
        {
            f[i][j] = invModFunc_(f[i][j]);
        }
    }
    return tmpf;
}


Foam::scalar Foam::lookupTable2D::lookup
(
    const scalar& x,
    const scalar& y
) const
{
    scalar fx, fy;
    label i, j;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    findYIndex_(modYFunc_(y), yMod_, j, fy);

    return
        invModFunc_
        (
            data_[i][j]*fx*fy
          + data_[i+1][j]*(1.0 - fx)*fy
          + data_[i][j+1]*fx*(1.0 - fy)
          + data_[i+1][j+1]*(1.0 - fx)*(1.0 - fy)
        );
}

Foam::scalar
Foam::lookupTable2D::reverseLookupY(const scalar& fin, const scalar& x) const
{
    scalar f(modFunc_(fin));
    scalar fx;
    label i;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    label j = boundj(f, data_, i);

    const scalar& mm(data_[i][j]);
    const scalar& pm(data_[i+1][j]);
    const scalar& mp(data_[i][j+1]);
    const scalar& pp(data_[i+1][j+1]);

    scalar fy =
        (f - fx*mp + fx*pp - pp)
       /(fx*mm - fx*mp - fx*pm + fx*pp + pm - pp);

    return invModYFunc_(getValue(j, 1.0 - fy, yMod_));
}


Foam::scalar
Foam::lookupTable2D::reverseLookupX(const scalar& fin, const scalar& y) const
{
    scalar f(modFunc_(fin));
    scalar fy;
    label j;
    findYIndex_(modYFunc_(y), yMod_, j, fy);
    label i = boundi(f, data_, j);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    scalar fx =
        (f - pm*fy - pp*(1.0 - fy))
       /(mm*fy + mp*(1.0 - fy) - (pm*fy + pp*(1.0 - fy)));

    return invModXFunc_(getValue(i, 1.0 - fx, xMod_));
}


Foam::scalar Foam::lookupTable2D::dFdX(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    findYIndex_(modYFunc_(y), yMod_, j, fy);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    return
        (
            invModFunc_(pm*fy + pp*(1.0 - fy))
          - invModFunc_(mm*fy + mp*(1.0 - fy))
        )/(x_[i+1] - x_[i]);
}

Foam::scalar Foam::lookupTable2D::dFdY(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    findYIndex_(modYFunc_(y), yMod_, j, fy);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    return
        (
            invModFunc_(mp*fx + pp*(1.0 - fx))
          - invModFunc_(mm*fx + pm*(1.0 - fx))
        )/(y_[j+1] - y_[j]);
}

Foam::scalar Foam::lookupTable2D::d2FdX2(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    findYIndex_(modYFunc_(y), yMod_, j, fy);

    if (i == 0)
    {
        i++;
    }

    scalar gmm(invModFunc_(data_[i-1][j]));
    scalar gm(invModFunc_(data_[i][j]));
    scalar gpm(invModFunc_(data_[i+1][j]));

    scalar gmp(invModFunc_(data_[i-1][j+1]));
    scalar gp(invModFunc_(data_[i][j+1]));
    scalar gpp(invModFunc_(data_[i+1][j+1]));

    const scalar& xm(x_[i-1]);
    const scalar& xi(x_[i]);
    const scalar& xp(x_[i+1]);

    scalar gPrimepm((gpm - gm)/(xp - xi));
    scalar gPrimemm((gm - gmm)/(xi - xm));
    scalar gPrimepp((gpp - gp)/(xp - xi));
    scalar gPrimemp((gp - gmp)/(xi - xm));

    return
        fy*(gPrimepm - gPrimemm)/(0.5*(xp - xm))
      + (1.0 - fy)*(gPrimepp - gPrimemp)/(0.5*(xp - xm));
}

Foam::scalar Foam::lookupTable2D::d2FdY2(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    findYIndex_(modYFunc_(y), yMod_, j, fy);

    if (j == 0)
    {
        j++;
    }

    scalar gmm(invModFunc_(data_[i][j-1]));
    scalar gm(invModFunc_(data_[i][j]));
    scalar gmp(invModFunc_(data_[i][j+1]));

    scalar gpm(invModFunc_(data_[i+1][j-1]));
    scalar gp(invModFunc_(data_[i+1][j]));
    scalar gpp(invModFunc_(data_[i+1][j+1]));

    const scalar& ym(y_[j-1]);
    const scalar& yi(y_[j]);
    const scalar& yp(y_[j+1]);

    scalar gPrimemp((gmp - gm)/(yp - yi));
    scalar gPrimemm((gm - gmm)/(yi - ym));
    scalar gPrimepp((gpp - gp)/(yp - yi));
    scalar gPrimepm((gp - gpm)/(yi - ym));

    return
        fx*(gPrimemp - gPrimemm)/(0.5*(yp - ym))
      + (1.0 - fx)*(gPrimepp - gPrimepm)/(0.5*(yp - ym));
}

Foam::scalar Foam::lookupTable2D::d2FdXdY(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findXIndex_(modXFunc_(x), xMod_, i, fx);
    findYIndex_(modYFunc_(y), yMod_, j, fy);

    scalar gmm(invModFunc_(data_[i][j]));
    scalar gmp(invModFunc_(data_[i][j+1]));
    scalar gpm(invModFunc_(data_[i+1][j]));
    scalar gpp(invModFunc_(data_[i+1][j+1]));

    const scalar& xm(x_[i]);
    const scalar& xp(x_[i+1]);

    const scalar& ym(y_[j]);
    const scalar& yp(y_[j+1]);

    return ((gpp - gmp)/(xp - xm) - (gpm - gmm)/(xp - xm))/(yp - ym);
}

// ************************************************************************* //
