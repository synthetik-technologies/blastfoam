/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "integrationRule.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::integrationRule::setTriMid
(
    const label start,
    const scalar w
)
{
    (*this)[start].set2(1.0/3.0, 1.0/3.0, w);
}


inline void Foam::integrationRule::setTri3
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar w
)
{
    UList<integrationPoint>& ip(*this);
    ip[start].set2(a, a, w);
    ip[start+1].set2(a, b, w);
    ip[start+2].set2(b, a, w);
}

inline void Foam::integrationRule::setTri3
(
    const label start,
    const scalar a,
    const scalar w
)
{
    setTri3(start, a, 1.0 - 2.0*a, w);
}

inline void Foam::integrationRule::setTri3b
(
    const label start,
    const scalar b,
    const scalar w
)
{
    setTri3(start, (1.0 - b)/2.0, b, w);
}

inline void Foam::integrationRule::setTri3R
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar w
)
{
    UList<integrationPoint>& ip(*this);
    ip[start].set2(a, b, w);
    ip[start+1].set2(c, a, w);
    ip[start+2].set2(b, c, w);
}

inline void Foam::integrationRule::setTri3R
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar w
)
{
    setTri3R(start, a, b, 1.0 - a - b, w);
}

inline void Foam::integrationRule::setTri6
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar w
)
{
    UList<integrationPoint>& ip(*this);
    ip[start].set2(a, b, w);
    ip[start+1].set2(b, a, w);
    ip[start+2].set2(a, c, w);
    ip[start+3].set2(c, a, w);
    ip[start+4].set2(b, c, w);
    ip[start+5].set2(c, b, w);
}

inline void Foam::integrationRule::setTri6
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar w
)
{
    setTri6(start, a, b, 1.0 - a - b, w);
}



inline void Foam::integrationRule::setTet3
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar w
)
{
    (*this)[start].set3(a, a, b, w);
    (*this)[start+1].set3(a, b, a, w);
    (*this)[start+2].set3(b, a, a, w);
}

inline void Foam::integrationRule::setTet6
(
    const label start,
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar w
)
{
    (*this)[start].set3(a, b, c, w);
    (*this)[start+1].set3(a, c, b, w);
    (*this)[start+2].set3(b, c, a, w);
    (*this)[start+3].set3(b, a, c, w);
    (*this)[start+4].set3(c, a, b, w);
    (*this)[start+5].set3(c, b, a, w);
}

inline void Foam::integrationRule::setTet6
(
    const label start,
    const scalar a,
    const scalar w
)
{
    const scalar b = 0.5 - a;
    setTet3(start, a, b, w);
    setTet3(start+3, b, a, w);
}

inline void Foam::integrationRule::setTetMid
(
    const label start,
    const scalar w
)
{
    (*this)[start].set3(0.25, 0.25, 0.25, w);
}

inline void Foam::integrationRule::setTet4
(
    const label start,
    const scalar a,
    const scalar w
)
{
    (*this)[start].set3(a, a, a, w);
    setTet3(start+1, a, 1.0 - 3.0*a, w);
}

inline void Foam::integrationRule::setTet4b
(
    const label start,
    const scalar b,
    const scalar w
)
{
    const scalar a = (1.0 - b)/3.0;
    (*this)[start].set3(a, a, a, w);
    setTet3(start+1, a, b, w);
}


inline void Foam::integrationRule::setTet12
(
    const label start,
    const scalar a,
    const scalar bc,
    const scalar w
)
{
    const scalar cb = 1.0 - 2.0*a - bc;
    setTet3(start, a, bc, w);
    setTet3(start+3, a, cb, w);
    setTet6(start+6, a, bc, cb, w);
}

inline void Foam::integrationRule::setTet12bc
(
    const label start,
    const scalar b,
    const scalar c,
    const scalar w
)
{
    const scalar a = (1.0 - b - c)/2.0;
    setTet3(start, a, b, w);
    setTet3(start+3, a, c, w);
    setTet6(start+6, a, b, c, w);
}

// ************************************************************************* //
