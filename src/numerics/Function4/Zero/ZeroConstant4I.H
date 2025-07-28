/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
01-11-2022 Synthetik Applied Technologies : Added Function4
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "ZeroConstant4.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function4s::ZeroConstant<Type>::value
(
    const scalar t,
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    return pTraits<Type>::zero;
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::Function4s::ZeroConstant<Type>::value
(
    const scalar t,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z
) const
{
    return tmp<Field<Type>>(new Field<Type>(x.size(), pTraits<Type>::zero));
}


// ************************************************************************* //
