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

#include "lookupTable1D.H"
#include "DynamicList.H"
#include "Field.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

void Foam::lookupTable1D::readTable()
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

    scalarField xTmp;
    scalarField xModTmp;
    scalarField fTmp;
    string line;
    while (is.good())
    {
        is.getLine(line);
        if (line[0] == '#')
        {
            continue;
        }

        IStringStream isLine(line);

        scalar xi(readScalar(isLine));

        if (is.good())
        {
            scalar fi(readScalar(isLine));
            xTmp.append(xi);
            xModTmp.append(modXFunc_(xi));
            fTmp.append(modFunc_(fi));
        }
    }
    xValues_ = xTmp;
    xModValues_ = xModTmp;
    data_ = fTmp;
}

void Foam::lookupTable1D::findIndex
(
    const scalar& x,
    label& i,
    scalar& f
) const
{
    if (x < xModValues_[0])
    {
        i = 0;
        f = x/xModValues_[0];
    }
    for (i = 1; i < xModValues_.size(); i++)
    {
        if (x < xModValues_[i])
        {
            i--;
            f = (x - xModValues_[i])/(xModValues_[i+1] - xModValues_[i]);
            return;
        }
    }
    i = xModValues_.size() - 2;
    f = 1.0;
    return;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupTable1D::lookupTable1D
(
    const fileName& file,
    const word& mod,
    const word& xMod
)
:
    file_(file),
    modFunc_(NULL),
    invModFunc_(NULL),
    modXFunc_(NULL),
    invModXFunc_(NULL)
{
    setMod(mod, modFunc_, invModFunc_);
    setMod(xMod, modXFunc_, invModXFunc_);
    readTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lookupTable1D::~lookupTable1D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::lookupTable1D::lookup(const scalar& x) const
{
    scalar fx;
    label i;
    findIndex(modXFunc_(x), i, fx);

    return invModFunc_(data_[i] + fx*(data_[i+1] - data_[i]));
}


Foam::scalar
Foam::lookupTable1D::reverseLookup(const scalar& fin) const
{
    scalar f(modFunc_(fin));
    if (f < data_[0])
    {
        return xValues_[0];
    }
    label i = 0;
    for (i = 1; i < data_.size(); i++)
    {
        if (f < data_[i])
        {
            i--;
            break;
        }
    }

    const scalar& fm(data_[i]);
    const scalar& fp(data_[i+1]);

    scalar fx = (f - fm)/(fp - fm);

    return invModXFunc_(xModValues_[i] + fx*(xValues_[i+1] - xValues_[i]));
}


Foam::scalar Foam::lookupTable1D::dFdX(const scalar& x) const
{
    scalar fx;
    label i;
    findIndex(modXFunc_(x), i, fx);

    scalar fm(data_[i]);
    scalar fp(data_[i+1]);

    return (invModFunc_(fp) - invModFunc_(fm))/(xValues_[i+1] - xValues_[i]);
}


Foam::scalar Foam::lookupTable1D::d2FdX2(const scalar& x) const
{
    label i;
    scalar fx;
    findIndex(modXFunc_(x), i, fx);

    if (i == 0)
    {
        i++;
    }

    scalar fm(invModFunc_(data_[i-1]));
    scalar fi(invModFunc_(data_[i]));
    scalar fp(invModFunc_(data_[i+1]));

    const scalar& xm(xValues_[i-1]);
    const scalar& xi(xValues_[i]);
    const scalar& xp(xValues_[i+1]);

    return
        ((fp - fi)/(xp - xi) - (fi - fm)/(xi - xm))/(xp - xm);
}

// ************************************************************************* //
