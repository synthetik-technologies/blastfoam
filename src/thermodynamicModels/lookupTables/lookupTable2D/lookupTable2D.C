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


void Foam::lookupTable2D::readTable()
{
    fileName fNameExpanded(file_);
    fNameExpanded.expand();

    // Open a stream and check it
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fNameExpanded));
    ISstream& is = isPtr();
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file" << file_ << nl
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
            data_[i][j] = readScalar(IStringStream(split[j])());
        }
        i++;
    }
}

void Foam::lookupTable2D::findIndex
(
    const scalar& xy,
    const scalar& xyMin,
    const scalar& dxy,
    const label nxy,
    label& IJ,
    scalar& f
) const
{
    scalar ij = (xy - xyMin)/dxy;
    if (ij < 0)
    {
//         if (debug)
//         {
//             WarningInFunction
//                 << name << " is out of bonds, limiting. "
//                 << "min(" << name << ") = "
//                 << invModVar(mod, getValue(0, xyMin, dxy))
//                 << ", " << name << " = " << xy << endl;
//         }
        IJ = 0;
        f = 1.0;
        return;
    }
    IJ = floor(ij);
    if (IJ >= nxy - 1)
    {
//         if ((ij >= (nxy - 1) && debug)
//         {
//             WarningInFunction
//                 << name << " is out of bonds, limiting. "
//                 << "max(" << name << " = "
//                 << invModVar(mod, getValue(ij, xyMin, dxy))
//                 << ", " << name << " = " << xy << endl;
//         }
        IJ = nxy - 2;
        f = 0.0;

        return;
    }

    f = 1.0 - (ij - scalar(IJ));
    return;
}

Foam::label Foam::lookupTable2D::bound
(
    const scalar& f,
    const label i,
    const scalar& yMin,
    const scalar& dy,
    const label ny,
    const bool ij
) const
{
    for (label j = 0; j < ny-1; j++)
    {
        if (ij)
        {
            if (f > data_[i][j] && f < data_[i+1][j])
            {
                return j;
            }
        }
        else
        {
            if (f > data_[j][i] && f < data_[j][i+1])
            {
                return j;
            }
        }
    }
    return ny - 2;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupTable2D::lookupTable2D
(
    const fileName& file,
    const word& mod,
    const word& xMod,
    const word& yMod,
    const label nx,
    const label ny,
    const scalar& xMin,
    const scalar& dx,
    const scalar& yMin,
    const scalar& dy
)
:
    file_(file),
    modFunc_(NULL),
    invModFunc_(NULL),
    modXFunc_(NULL),
    invModXFunc_(NULL),
    modYFunc_(NULL),
    invModYFunc_(NULL),
    nx_(nx),
    ny_(ny),
    xMin_(xMin),
    dx_(dx),
    yMin_(yMin),
    dy_(dy),
    data_(nx_, scalarField(ny_, 0.0)),
    x_(nx_, 0.0),
    y_(ny_, 0.0)
{
    if (mod == "log10")
    {
        modFunc_ = &log10S;
        invModFunc_ = &pow10S;
    }
    else if (mod == "ln")
    {
        modFunc_ = &lnS;
        invModFunc_ = &expS;
    }
    else if (mod == "exp")
    {
        modFunc_ = &expS;
        invModFunc_ = &lnS;
    }
    else
    {
        modFunc_ = &noneS;
        invModFunc_ = &noneS;
    }

    if (xMod == "log10")
    {
        modXFunc_ = &log10S;
        invModXFunc_ = &pow10S;
    }
    else if (xMod == "ln")
    {
        modXFunc_ = &lnS;
        invModXFunc_ = &expS;
    }
    else if (xMod ==  "exp")
    {
        modXFunc_ = &expS;
        invModXFunc_ = &lnS;
    }
    else
    {
        modXFunc_ = &noneS;
        invModXFunc_ = &noneS;
    }

    if (yMod == "log10")
    {
        modYFunc_ = &log10S;
        invModYFunc_ = &pow10S;
    }
    else if (yMod == "ln")
    {
        modYFunc_ = &lnS;
        invModYFunc_ = &expS;
    }
    else if (yMod ==  "exp")
    {
        modYFunc_ = &expS;
        invModYFunc_ = &lnS;
    }
    else
    {
        modYFunc_ = &noneS;
        invModYFunc_ = &noneS;
    }

    readTable();

    forAll(x_, i)
    {
        x_[i] = invModXFunc_(getValue(i, xMin_, dx_));
    }
    forAll(y_, i)
    {
        y_[i] = invModYFunc_(getValue(i, yMin_, dy_));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lookupTable2D::~lookupTable2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::lookupTable2D::lookup
(
    const scalar& x,
    const scalar& y
) const
{
    scalar fx, fy;
    label i, j;
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);

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
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    label j = bound(f, i, yMin_, dy_, ny_, true);

    const scalar& mm(data_[i][j]);
    const scalar& pm(data_[i+1][j]);
    const scalar& mp(data_[i][j+1]);
    const scalar& pp(data_[i+1][j+1]);

    scalar fy =
        (f - fx*mp + fx*pp - pp)
       /(fx*mm - fx*mp - fx*pm + fx*pp + pm - pp);

    return invModYFunc_(getValue(j + 1.0 - fy, yMin_, dy_));
}


Foam::scalar
Foam::lookupTable2D::reverseLookupX(const scalar& fin, const scalar& y) const
{
    scalar f(modFunc_(fin));
    scalar fy;
    label j;
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);
    label i = bound(f, j, xMin_, dx_, nx_, false);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    scalar fx =
        (f - pm*fy - pp*(1.0 - fy))
       /(mm*fy + mp*(1.0 - fy) - (pm*fy + pp*(1.0 - fy)));

    return invModXFunc_(getValue(i + 1.0 - fx, xMin_, dx_));
}


Foam::scalar Foam::lookupTable2D::dFdX(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);

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
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);

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
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);

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
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);

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
    findIndex(modXFunc_(x), xMin_, dx_, nx_, i, fx);
    findIndex(modYFunc_(y), yMin_, dy_, ny_, j, fy);

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
