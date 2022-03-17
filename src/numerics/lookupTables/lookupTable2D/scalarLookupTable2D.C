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

#include "scalarLookupTable2D.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

Foam::labelList Foam::scalarLookupTable2D::boundi
(
    const scalar f
) const
{
    const label i = ij_.x();
    const label j = ij_.y();
    if (f < data_[0][j])
    {
        return labelList(1, 0);
    }

    DynamicList<label> Is(data_.size());
    for (label I = 0; I < data_.size() - 1; I++)
    {
        if
        (
            f > data_[i][j]
         && f < data_[i][j+1]
         && f > data_[i+1][j]
         && f < data_[i+1][j+1]
        )
        {
            Is.append(I);
        }
    }
    if (!Is.size())
    {
        return labelList(1, data_.size() - 2);
    }
    return move(Is);
}


Foam::labelList Foam::scalarLookupTable2D::boundj
(
    const scalar f
) const
{
    const label i = ij_.x();
    const label j = ij_.y();
    if (data_[i][0] > f)
    {
        return labelList(1, 0);
    }

    DynamicList<label> Js(data_[i].size());
    for (label J = 0; J < data_[i].size() - 1; J++)
    {
        if
        (
            f > data_[i][j]
         && f < data_[i+1][j]
         && f > data_[i][j+1]
         && f < data_[i+1][j+1]
        )
        {
            Js.append(j);
        }
    }
    if (!Js.size())
    {
        return labelList(1, data_[i].size() - 2);
    }
    return move(Js);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarLookupTable2D::scalarLookupTable2D()
{}


Foam::scalarLookupTable2D::scalarLookupTable2D
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name,
    const bool canRead
)
{
    read(dict, xName, yName, name, canRead);
}


Foam::scalarLookupTable2D::scalarLookupTable2D
(
    const Field<scalar>& x,
    const Field<scalar>& y,
    const Field<Field<scalar>>& data,
    const word& modXType,
    const word& modYType,
    const word& modType,
    const word& interpolationScheme,
    const bool isReal
)
{
    set
    (
        x,
        y,
        data,
        modXType,
        modYType,
        modType,
        interpolationScheme,
        isReal
    );
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

    ij_.x() = findXIndex_(x, xValues());
    labelList Js(boundj(f));

    const label i = ij_.x();
    label& j = ij_.y();
    scalar fx = linearWeight(modXFunc_(x), xModValues_[i], xModValues_[i+1]);
    scalar fy = 1.0;

    if (Js.size() == 1)
    {
        j = Js[0];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);
        fy =
            (f + fx*(mm  - pm) - mm)
           /(fx*(mm - pm - mp + pp) - mm + mp);
        return invModYFunc_(getValue(j, fy, yModValues_));
    }

    //- If multiple indicies meet criteria, check for closest
    Field<scalar> yTrys(Js.size());
    Field<scalar> errors(Js.size(), great);
    forAll(Js, J)
    {
        j = Js[J];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);
        fy =
            (f + fx*(mm  - pm) - mm)
           /(fx*(mm - pm - mp + pp) - mm + mp);
        yTrys[J] = invModYFunc_(getValue(j, fy, yModValues_));
        errors[J] = mag(fin - lookup(x, yTrys[j]));
    }
    j = findMin(errors);
    return yTrys[j];
}


Foam::scalar
Foam::scalarLookupTable2D::reverseLookupX
(
    const scalar fin,
    const scalar y
) const
{
    scalar f(modFunc_(fin));
    ij_.y() = findYIndex_(y, yValues());
    labelList Is(boundi(f));

    const label j = ij_.y();
    label& i = ij_.x();
    scalar fy = linearWeight(modYFunc_(y), yModValues_[j], yModValues_[j+1]);
    scalar fx = 1.0;

    if (Is.size() == 1)
    {
        i = Is[0];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);

        fx =
            (f + fy*(mm - mp) - mm)
           /(fy*(mm - mp - pm + pp) - mm + pm);

        return invModXFunc_(getValue(i, fx, xModValues_));
    }

    //- If multiple indicies meet criteria, check for closest
    Field<scalar> xTrys(Is.size());
    Field<scalar> errors(Is.size(), great);
    forAll(Is, I)
    {
        i = Is[I];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);
        fx =
            (f + fy*(mm - mp) - mm)
           /(fy*(mm - mp - pm + pp) - mm + pm);
        xTrys[I] = invModXFunc_(getValue(i, fx, xModValues_));
        errors[I] = mag(fin - lookup(xTrys[I], y));
    }
    i = findMin(errors);
    return xTrys[i];
}

// ************************************************************************* //
