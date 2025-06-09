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

#include "integrationPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::integrationPoint::integrationPoint()
:
    vector(0.0, 0.0, 0.0),
    w_(0.0)
{}


inline Foam::integrationPoint::integrationPoint
(
    const vector& pt,
    const scalar w
)
:
    vector(pt),
    w_(w)
{}


inline Foam::integrationPoint::integrationPoint
(
    const scalar x,
    const scalar y,
    const scalar z,
    const scalar w
)
:
    vector(x, y, z),
    w_(w)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::integrationPoint::~integrationPoint()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::integrationPoint::w() const
{
    return w_;
}


inline Foam::scalar& Foam::integrationPoint::w()
{
    return w_;
}


inline void Foam::integrationPoint::set1(const scalar x)
{
    v_[0] = x;
}


inline void Foam::integrationPoint::set1(const scalar x, const scalar w)
{
    v_[0] = x;
    w_ = w;
}


inline void Foam::integrationPoint::set2
(
    const scalar x,
    const scalar y
)
{
    v_[0] = x;
    v_[1] = y;
}


inline void Foam::integrationPoint::set2
(
    const scalar x,
    const scalar y,
    const scalar w
)
{
    v_[0] = x;
    v_[1] = y;
    w_ = w;
}


inline void Foam::integrationPoint::set3
(
    const scalar x,
    const scalar y,
    const scalar z
)
{
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
}


inline void Foam::integrationPoint::set3
(
    const scalar x,
    const scalar y,
    const scalar z,
    const scalar w
)
{
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
    w_ = w;
}

// ************************************************************************* //
