/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    HormannAgathos

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Implements a == b within a tolerance define by epsilon_.
// Need to work if epsilon_ == 0.0
inline bool Foam::HormannAgathos::equalWithTol
(
    const scalar& a,
    const scalar& b
) const
{
    return (a <= b + epsilon_) && (a >= b - epsilon_);
}


// Implements a > b within a tolerance define by epsilon_.
// Need to work if epsilon_ == 0.0
inline bool Foam::HormannAgathos::greaterWithTol
(
    const scalar& a,
    const scalar& b
) const
{
    return a > (b + epsilon_);
}


// Implements a < b within a tolerance define by epsilon_.
// Need to work if epsilon_ == 0.0
inline bool Foam::HormannAgathos::smallerWithTol
(
    const scalar& a,
    const scalar& b
) const
{
    return a < (b - epsilon_);
}


// Implements a >= b within a tolerance define by epsilon_.
// Need to work if epsilon_ == 0.0
inline bool Foam::HormannAgathos::greaterOrEqualWithTol
(
    const scalar& a,
    const scalar& b
) const
{
    return a >= (b - epsilon_);
}


// Implements a =< b within a tolerance define by epsilon_.
// Need to work if epsilon_ == 0.0
inline bool Foam::HormannAgathos::smallerOrEqualWithTol
(
    const scalar& a,
    const scalar& b
) const
{
    return a <= (b + epsilon_);
}


// ************************************************************************* //
