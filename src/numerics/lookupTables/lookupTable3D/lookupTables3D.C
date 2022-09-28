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

#include "lookupTables3D.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y,
    const scalar z
) const
{
    return solver(0, fin, y, z).solveUni();
}


template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x,
    const scalar z
) const
{
    return solver(1, fin, x, z).solveUni();
}


template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupZ
(
    const scalar& fin,
    const scalar x,
    const scalar y
) const
{

    return solver(2, fin, x, y).solveUni();
}

// ************************************************************************* //
