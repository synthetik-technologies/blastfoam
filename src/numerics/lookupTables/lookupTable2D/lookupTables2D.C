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

#include "lookupTables2D.H"

// * * * * * * * * * * * * * * Scalar Functions * * * * * * * * * * * * * * //

template<>
Foam::labelList Foam::lookupTable2D<Foam::scalar>::boundi
(
    const scalar& f
) const
{
    label& i = ij_.x();
    const label j = ij_.y();
    if (f < data_[0][j])
    {
        return labelList(1, 0);
    }

    DynamicList<label> Is(data_.size());
    for (i = 0; i < data_.size(); i++)
    {
        if (f > data_[i][j] && f < data_[i][j+1])
        {
            Is.append(i);
        }
    }
    if (!Is.size())
    {
        return labelList(1, data_.size() - 2);
    }
    if (Is.size() == 1)
    {
        return move(Is);
    }
    DynamicList<label> newIs(Is.size());
    for (i = 0; i < Is.size() - 1; i++)
    {
        if (Is[i]+1 == Is[i+1])
        {
            newIs.append(Is[i]);
        }
    }
    return move(newIs);
}


template<>
Foam::labelList Foam::lookupTable2D<Foam::scalar>::boundj
(
    const scalar& f
) const
{
    const label i = ij_.x();
    label& j = ij_.y();
    if (data_[i][0] > f)
    {
        return labelList(1, 0);
    }

    DynamicList<label> Js(data_[i].size());
    for (j = 0; j < data_[i].size(); j++)
    {
        if (f > data_[i][j] && f < data_[i+1][j])
        {
            Js.append(j);
        }
    }
    if (!Js.size())
    {
        return labelList(1, data_[i].size() - 2);
    }
    if (Js.size() == 1)
    {
        return move(Js);
    }
    DynamicList<label> newJs(Js.size());
    for (j = 0; j < Js.size() - 1; j++)
    {
        if (Js[j]+1 == Js[j+1])
        {
            newJs.append(Js[j]);
        }
    }
    return move(newJs);
}


template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y
) const
{
    scalar f(mod_()(fin));
    ij_.y() = yIndexing_->findIndex(modY_()(y));
    labelList Is(boundi(f));

    const label j = ij_.y();
    label& i = ij_.x();
    scalar fy = interpolationWeight1D::linearWeight
    (
        modY_()(y),
        yModValues_[j],
        yModValues_[j+1]
    );
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

        return modX_->inv(getValue(i, fx, xModValues_));
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
        xTrys[I] = modX_->inv(getValue(i, fx, xModValues_));
        errors[I] = mag(fin - lookup(xTrys[I], y));
    }
    i = findMin(errors);
    return xTrys[i];
}


template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x
) const
{
    scalar f(mod_()(fin));

    ij_.x() = xIndexing_->findIndex(modX_()(x));
    labelList Js(boundj(f));

    const label i = ij_.x();
    label& j = ij_.y();
    scalar fx = interpolationWeight1D::linearWeight
    (
        modX_()(x),
        xModValues_[i],
        xModValues_[i+1]
    );
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
        return modY_->inv(getValue(j, fy, yModValues_));
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
        yTrys[J] = modY_->inv(getValue(j, fy, yModValues_));
        errors[J] = mag(fin - lookup(x, yTrys[J]));
    }
    j = findMin(errors);
    return yTrys[j];
}


// ************************************************************************* //
