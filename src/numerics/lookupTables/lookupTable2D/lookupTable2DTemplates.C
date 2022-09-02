/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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

#include "lookupTable2D.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<template<class> class ListType, class fType>
fType Foam::lookupTable2D<Type>::interpolate
(
    const List<ListType<fType>>& fs
) const
{
    fType modf = weights_[0]*fs[indices_[0].x()][indices_[0].y()];
    for (label i = 1; i < indices_.size(); i++)
    {
        modf += weights_[i]*fs[indices_[i].x()][indices_[i].y()];
    }
    return modf;
}


template<class Type>
template<template<class> class ListType, class fType>
fType Foam::lookupTable2D<Type>::interpolate
(
    const scalar x,
    const scalar y,
    const List<ListType<fType>>& fs
) const
{
    update(x, y);

    return interpolate(fs);
}

// ************************************************************************* //
