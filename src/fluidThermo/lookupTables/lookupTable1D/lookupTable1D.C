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
#include "IFstream.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

void Foam::lookupTable1D::readTable()
{
    IFstream is(file_);
    string line;
    bool good = is.good();
    if (!good)
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file" << file_ << nl
            << exit(FatalIOError);
    }
    while (good)
    {
        scalar xTmp(readScalar(is));
        good = is.good();
        if (good)
        {
            x_.append(modXFunc_(xTmp));
            realX_.append(xTmp);
            data_.append(modFunc_(readScalar(is)));
        }
    }
}

void Foam::lookupTable1D::findIndex
(
    const scalar& x,
    label& I,
    scalar& f
) const
{
    for (label i = 0; i < x_.size(); i++)
    {
        if (x_[i] > x)
        {
            I = max(i-1, 0);
            f = (x - x_[I])/(x_[I+1] - x_[I]);
            return;
        }
    }
    I = x_.size() - 2;
    f = 0.0;
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

    return invModFunc_(data_[i]*fx + data_[i+1]*(1.0 - fx));
}


Foam::scalar
Foam::lookupTable1D::reverseLookup(const scalar& fin) const
{
    scalar f(modFunc_(fin));
    label i;
    for (i = 0; i <data_.size(); i++)
    {
        if (data_[i] < f)
        {
            break;
        }
    }

    scalar fx = (data_[i+1] - fin)/(x_[i+1] - x_[i]);

    return invModXFunc_(x_[i] + (1.0 - fx)*(x_[i+1] - x_[i]));
}


Foam::scalar Foam::lookupTable1D::dFdX(const scalar& x) const
{
    scalar fx;
    label i;
    findIndex(modXFunc_(x), i, fx);

    return (invModFunc_(data_[i+1]) - invModFunc_(data_[i]))/(realX_[i+1] - realX_[i]);
}

Foam::scalar Foam::lookupTable1D::d2FdX2(const scalar& x) const
{
    scalar fx;
    label i;
    findIndex(modXFunc_(x), i, fx);

    if (i == 0)
    {
        i++;
    }

    scalar gm(invModFunc_(data_[i-1]));
    scalar g(invModFunc_(data_[i]));
    scalar gp(invModFunc_(data_[i+1]));

    const scalar& xm(realX_[i-1]);
    const scalar& xi(realX_[i]);
    const scalar& xp(realX_[i+1]);

    return (gp - 2.0*g + gm)/((xp - xi)*(xi - xm));
}

// ************************************************************************* //
