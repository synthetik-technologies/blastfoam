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

#include "scalarLookupTable2D.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

Foam::labelList Foam::scalarLookupTable2D::boundi
(
    const scalar f
) const
{
    if (f < data_[0][j_])
    {
        return labelList(1, 0);
    }

    labelList I(data_.size());
    label nFound = 0;
    for (i_ = 0; i_ < data_.size() - 1; i_++)
    {
        if
        (
            f > data_[i_][j_]
         && f < data_[i_][j_+1]
         && f > data_[i_+1][j_]
         && f < data_[i_+1][j_+1]
        )
        {
            I[nFound++] = i_;
        }
    }
    if (!nFound)
    {
        return labelList(1, data_.size() - 2);
    }
    I.resize(nFound);
    return I;
}


Foam::labelList Foam::scalarLookupTable2D::boundj
(
    const scalar f
) const
{
    if (data_[i_][0] > f)
    {
        return labelList(1, 0);
    }

    labelList J(data_[i_].size());
    label nFound = 0;
    for (j_ = 0; j_ < data_[i_].size() - 1; j_++)
    {
        if
        (
            f > data_[i_][j_]
         && f < data_[i_+1][j_]
         && f > data_[i_][j_+1]
         && f < data_[i_+1][j_+1]
        )
        {
            J[nFound++] = j_;
        }
    }
    if (!nFound)
    {
        return labelList(1, data_[i_].size() - 2);
    }
    J.resize(nFound);
    return J;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarLookupTable2D::scalarLookupTable2D()
{}


Foam::scalarLookupTable2D::scalarLookupTable2D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name
)
{
    read(dict, xName, yName, name);
}


Foam::scalarLookupTable2D::scalarLookupTable2D
(
    const scalarField& x,
    const scalarField& y,
    const Field<scalarField>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const bool isReal
)
{
    set(x, y, data, modXType, modYType, modType, isReal);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarLookupTable2D::~scalarLookupTable2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::scalarLookupTable2D::reverseLookupY
(
    const scalar fin,
    const scalar x
) const
{
    scalar f(modFunc_(fin));
    updateX(x);
    labelList Js(boundj(f));

    if (Js.size() == 1)
    {
        j_ = Js[0];
        const scalar mm(data_[i_][j_]);
        const scalar pm(data_[i_+1][j_]);
        const scalar mp(data_[i_][j_+1]);
        const scalar pp(data_[i_+1][j_+1]);
        fy_ =
            (f - fx_*mp + fx_*pp - pp)
           /(fx_*mm - fx_*mp - fx_*pm + fx_*pp + pm - pp);
        return invModYFunc_(getValue(j_, 1.0 - fy_, yModValues_));
    }

    //- If multiple indicies meet criteria, check for closest
    scalarList yTrys(Js.size());
    scalarList errors(Js.size(), great);
    forAll(Js, J)
    {
        j_ = Js[J];
        const scalar mm(data_[i_][j_]);
        const scalar pm(data_[i_+1][j_]);
        const scalar mp(data_[i_][j_+1]);
        const scalar pp(data_[i_+1][j_+1]);
        fy_ =
            (f - fx_*mp + fx_*pp - pp)
           /(fx_*mm - fx_*mp - fx_*pm + fx_*pp + pm - pp);
        yTrys[J] = invModYFunc_(getValue(j_, 1.0 - fy_, yModValues_));
        errors[J] = mag(fin - lookup(x, yTrys[j_]));
    }
    j_ = findMin(errors);
    return yTrys[j_];
}


Foam::scalar
Foam::scalarLookupTable2D::reverseLookupX
(
    const scalar fin,
    const scalar y
) const
{
    scalar f(modFunc_(fin));
    updateY(y);
    labelList Is(boundi(f));

    if (Is.size() == 1)
    {
        i_ = Is[0];
        const scalar mm(data_[i_][j_]);
        const scalar pm(data_[i_+1][j_]);
        const scalar mp(data_[i_][j_+1]);
        const scalar pp(data_[i_+1][j_+1]);
        fx_ =
            (f - pm*fy_ - pp*(1.0 - fy_))
           /(mm*fy_ + mp*(1.0 - fy_) - (pm*fy_ + pp*(1.0 - fy_)));
        return invModXFunc_(getValue(i_, 1.0 - fx_, xModValues_));
    }

    //- If multiple indicies meet criteria, check for closest
    scalarList xTrys(Is.size());
    scalarList errors(Is.size(), great);
    forAll(Is, I)
    {
        i_ = Is[I];
        const scalar mm(data_[i_][j_]);
        const scalar pm(data_[i_+1][j_]);
        const scalar mp(data_[i_][j_+1]);
        const scalar pp(data_[i_+1][j_+1]);
        fx_ =
            (f - pm*fy_ - pp*(1.0 - fy_))
           /(mm*fy_ + mp*(1.0 - fy_) - (pm*fy_ + pp*(1.0 - fy_)));
        xTrys[I] = invModXFunc_(getValue(i_, 1.0 - fx_, xModValues_));
        errors[I] = mag(fin - lookup(xTrys[I], y));
    }
    i_ = findMin(errors);
    return xTrys[i_];
}

// ************************************************************************* //
