/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "lookupTable3D.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template
<
    class fType,
    template<class> class ListType1,
    template<class> class ListType2
>
fType Foam::lookupTable3D<Type>::interpolate
(
    const List<ListType1<ListType2<fType>>>& fs
) const
{
    fType modf =
        weights_[0]
       *fs[indices_[0].x()][indices_[0].y()][indices_[0].z()];
    for (label i = 1; i < indices_.size(); i++)
    {
        modf +=
            weights_[i]
           *fs[indices_[i].x()][indices_[i].y()][indices_[i].z()];
    }
    return modf;
}


template<class Type>
template
<
    class fType,
    template<class> class ListType1,
    template<class> class ListType2
>
fType Foam::lookupTable3D<Type>::interpolate
(
    const scalar x,
    const scalar y,
    const scalar z,
    const List<ListType1<ListType2<fType>>>& fs
) const
{
    update(x, y, z);

    return interpolate(fs);
}

// ************************************************************************* //
