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

#include "lookupTable.H"
#include "DynamicList.H"
#include "Field.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

Foam::lookupTable::modType
Foam::lookupTable::getModType(const word& type) const
{
    if (type == "log10")
    {
        return modType::log10;
    }
    else if (type == "ln")
    {
        return modType::ln;
    }
    else if (type == "exp")
    {
        return modType::exp;
    }
    return modType::none;
}


Foam::scalar Foam::lookupTable::readValue(const List<string>& split) const
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


void Foam::lookupTable::readTable()
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

void Foam::lookupTable::findIndex
(
    const word& name,
    const scalar& xy,
    const scalar& xyMin,
    const scalar& dxy,
    const label nxy,
    const modType& mod,
    label& IJ,
    scalar& f
) const
{
    scalar ij = (modVar(mod, xy) - xyMin)/dxy;
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

Foam::label Foam::lookupTable::bound
(
    const scalar& f,
    const label i,
    const scalar& yMin,
    const scalar& dy,
    const label ny,
    const bool ij
) const
{
    for (label j = 0; j < ny; j++)
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

Foam::scalar
Foam::lookupTable::modVar(const modType& type, const scalar& xy) const
{
    if (type == log10)
    {
        return Foam::log10(max(xy, small));
    }
    else if (type == ln)
    {
        return Foam::log(max(xy, small));
    }
    else if (type == exp)
    {
        return Foam::exp(xy);
    }

    return xy;
}

Foam::scalar
Foam::lookupTable::invModVar(const modType& type, const scalar& xy) const
{
    if (type == log10)
    {
        return Foam::pow(10.0, xy);
    }
    else if (type == ln)
    {
        return Foam::exp(xy);
    }
    else if (type == exp)
    {
        return Foam::log(max(xy, small));
    }

    return xy;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupTable::lookupTable
(
    const word& fileName,
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
    file_(fileName),
    mod_(getModType(mod)),
    xMod_(getModType(xMod)),
    yMod_(getModType(yMod)),
    nx_(nx),
    ny_(ny),
    xMin_(xMin),
    dx_(dx),
    yMin_(yMin),
    dy_(dy),
    data_(nx_, scalarField(ny_, 0.0))
{
    readTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lookupTable::~lookupTable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::lookupTable::lookup
(
    const scalar& x,
    const scalar& y
) const
{
    scalar fx, fy;
    label i, j;
    findIndex("x", x, xMin_, dx_, nx_, xMod_, i, fx);
    findIndex("y", y, yMin_, dy_, ny_, yMod_, j, fy);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    return
        invModVar
        (
            mod_,
            mm*fx*fy + pm*(1.0 - fx)*fy + mp*fx*(1.0 - fy) + pp*(1.0 - fx)*(1.0 - fy)
        );
}

Foam::scalar
Foam::lookupTable::reverseLookupY(const scalar& fin, const scalar& x) const
{
    scalar f(modVar(mod_, fin));
    scalar fx;
    label i;
    findIndex("x", max(x, small), xMin_, dx_, nx_, xMod_, i, fx);
    label j = bound(f, i, yMin_, dy_, ny_, true);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    scalar fy =
        (f - fx*mp + fx*pp - pp)
       /(fx*mm - fx*mp - fx*pm + fx*pp + pm - pp);

    return invModVar(yMod_, getValue(j + 1.0 - fy, yMin_, dy_));
}


Foam::scalar
Foam::lookupTable::reverseLookupX(const scalar& fin, const scalar& y) const
{
    scalar f(modVar(mod_, fin));
    scalar fy;
    label j;
    findIndex("y", max(y, small), yMin_, dy_, ny_, yMod_, j, fy);
    label i = bound(f, j, xMin_, dx_, nx_, false);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    scalar fx =
        (f - pm*fy - pp*(1.0 - fy))
       /(mm*fy + mp*(1.0 - fy) - (pm*fy + pp*(1.0 - fy)));

    return invModVar(xMod_, getValue(i + 1.0 - fx, xMin_, dx_));
}


Foam::scalar Foam::lookupTable::dFdX(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findIndex("x", max(x, small), xMin_, dx_, nx_, xMod_, i, fx);
    findIndex("y", max(y, small), yMin_, dy_, ny_, yMod_, j, fy);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    return
        (
            invModVar(mod_, pm*fy + pp*(1.0 - fy))
          - invModVar(mod_, mm*fy + mp*(1.0 - fy))
        )
       /(
           invModVar(xMod_, getValue(i+1, xMin_, dx_))
         - invModVar(xMod_, getValue(i, xMin_, dx_))
        );
}

Foam::scalar Foam::lookupTable::dFdY(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findIndex("x", max(x, small), xMin_, dx_, nx_, xMod_, i, fx);
    findIndex("y", max(y, small), yMin_, dy_, ny_, yMod_, j, fy);

    scalar mm(data_[i][j]);
    scalar pm(data_[i+1][j]);
    scalar mp(data_[i][j+1]);
    scalar pp(data_[i+1][j+1]);

    return
        (
            invModVar(mod_, mp*fx + pp*(1.0 - fx))
          - invModVar(mod_, mm*fx + pm*(1.0 - fx))
        )
       /(
           invModVar(yMod_, getValue(j+1, yMin_, dy_))
         - invModVar(yMod_, getValue(j, yMin_, dy_))
        );
}

Foam::scalar Foam::lookupTable::d2FdX2(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findIndex("x", max(x, small), xMin_, dx_, nx_, xMod_, i, fx);
    findIndex("y", max(y, small), yMin_, dy_, ny_, yMod_, j, fy);

    if (i == 0)
    {
        i++;
    }

    scalar gmm(invModVar(mod_, data_[i-1][j]));
    scalar gm(invModVar(mod_, data_[i][j]));
    scalar gpm(invModVar(mod_, data_[i+1][j]));

    scalar gmp(invModVar(mod_, data_[i-1][j+1]));
    scalar gp(invModVar(mod_, data_[i][j+1]));
    scalar gpp(invModVar(mod_, data_[i+1][j+1]));

    scalar xm(invModVar(xMod_, getValue(i-1, xMin_, dx_)));
    scalar xi(invModVar(xMod_, getValue(i, xMin_, dx_)));
    scalar xp(invModVar(xMod_, getValue(i+1, xMin_, dx_)));

    scalar gPrimepm((gpm - gm)/(xp - xi));
    scalar gPrimemm((gm - gmm)/(xi - xm));
    scalar gPrimepp((gpp - gp)/(xp - xi));
    scalar gPrimemp((gp - gmp)/(xi - xm));

    return
        fy*(gPrimepm - gPrimemm)/(0.5*(xp - xm))
      + (1.0 - fy)*(gPrimepp - gPrimemp)/(0.5*(xp - xm));
}

Foam::scalar Foam::lookupTable::d2FdXdY(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findIndex("x", max(x, small), xMin_, dx_, nx_, xMod_, i, fx);
    findIndex("y", max(y, small), yMin_, dy_, ny_, yMod_, j, fy);

    if (j == 0)
    {
        j++;
    }

    scalar gmm(invModVar(mod_, data_[i][j-1]));
    scalar gm(invModVar(mod_, data_[i][j]));
    scalar gmp(invModVar(mod_, data_[i][j+1]));

    scalar gpm(invModVar(mod_, data_[i+1][j-1]));
    scalar gp(invModVar(mod_, data_[i+1][j]));
    scalar gpp(invModVar(mod_, data_[i+1][j+1]));

    scalar ym(invModVar(yMod_, getValue(j-1, yMin_, dy_)));
    scalar yi(invModVar(yMod_, getValue(j, yMin_, dy_)));
    scalar yp(invModVar(yMod_, getValue(j+1, yMin_, dy_)));

    scalar gPrimemp((gmp - gm)/(yp - yi));
    scalar gPrimemm((gm - gmm)/(yi - ym));
    scalar gPrimepp((gpp - gp)/(yp - yi));
    scalar gPrimepm((gp - gpm)/(yi - ym));

    return
        fx*(gPrimemp - gPrimemm)/(0.5*(yp - ym))
      + (1.0 - fx)*(gPrimepp - gPrimepm)/(0.5*(yp - ym));
}

Foam::scalar Foam::lookupTable::d2FdY2(const scalar& x, const scalar& y) const
{
    scalar fx, fy;
    label i, j;
    findIndex("x", max(x, small), xMin_, dx_, nx_, xMod_, i, fx);
    findIndex("y", max(y, small), yMin_, dy_, ny_, yMod_, j, fy);

    if (j == 0)
    {
        j++;
    }
    if (i == 0)
    {
        i++;
    }

    scalar gmm(invModVar(mod_, data_[i][j]));
    scalar gmp(invModVar(mod_, data_[i][j+1]));
    scalar gpm(invModVar(mod_, data_[i+1][j]));
    scalar gpp(invModVar(mod_, data_[i+1][j+1]));

    scalar xm(invModVar(xMod_, getValue(i, xMin_, dx_)));
    scalar xp(invModVar(xMod_, getValue(i+1, xMin_, dx_)));

    scalar ym(invModVar(yMod_, getValue(j, yMin_, dy_)));
    scalar yp(invModVar(yMod_, getValue(j+1, yMin_, dy_)));

    return ((gpp - gmp)/(xp - xm) - (gpm - gmm)/(xp - xm))/(yp - ym);
}

// ************************************************************************* //
