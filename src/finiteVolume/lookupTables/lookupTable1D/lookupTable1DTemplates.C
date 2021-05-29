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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::lookupTable1D::interpolate
(
    const Type& fm,
    const Type& fp
) const
{
    return f_ == 0 ? fm : (f_ == 1 ? fp : (1.0 - f_)*fm + f_*fp);
}


template<template<class> class ListType, class Type>
Type Foam::lookupTable1D::interpolate
(
    const scalar x,
    const ListType<Type>& fs
) const
{
    update(x);
    return
        f_ == 0 ? fs[index_]
      : (
            f_ == 1
          ? fs[index_ + 1]
          : (1.0 - f_)*fs[index_] + f_*fs[index_ + 1]
        );
}

// ************************************************************************* //
