/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::objectMotionState::centreOfRotation() const
{
    return centreOfRotation_;
}


inline const Foam::tensor& Foam::objectMotionState::Q() const
{
    return Q_;
}


inline const Foam::vector& Foam::objectMotionState::v() const
{
    return v_;
}


inline const Foam::vector& Foam::objectMotionState::a() const
{
    return a_;
}


inline const Foam::vector& Foam::objectMotionState::pi() const
{
    return pi_;
}


inline const Foam::vector& Foam::objectMotionState::tau() const
{
    return tau_;
}


inline Foam::point& Foam::objectMotionState::centreOfRotation()
{
    return centreOfRotation_;
}


inline Foam::tensor& Foam::objectMotionState::Q()
{
    return Q_;
}


inline Foam::vector& Foam::objectMotionState::v()
{
    return v_;
}


inline Foam::vector& Foam::objectMotionState::a()
{
    return a_;
}


inline Foam::vector& Foam::objectMotionState::pi()
{
    return pi_;
}


inline Foam::vector& Foam::objectMotionState::tau()
{
    return tau_;
}


// ************************************************************************* //
