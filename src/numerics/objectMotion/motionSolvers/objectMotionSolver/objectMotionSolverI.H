/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

inline Foam::point& Foam::objectMotionSolver::centreOfRotation()
{
    return motionState_.centreOfRotation();
}

inline Foam::tensor& Foam::objectMotionSolver::Q()
{
    return motionState_.Q();
}

inline Foam::vector& Foam::objectMotionSolver::v()
{
    return motionState_.v();
}

inline Foam::vector& Foam::objectMotionSolver::a()
{
    return motionState_.a();
}

inline Foam::vector& Foam::objectMotionSolver::pi()
{
    return motionState_.pi();
}

inline Foam::vector& Foam::objectMotionSolver::tau()
{
    return motionState_.tau();
}


inline const Foam::point& Foam::objectMotionSolver::centreOfRotation0() const
{
    return motionState0_.centreOfRotation();
}

inline const Foam::tensor& Foam::objectMotionSolver::Q0() const
{
    return motionState0_.Q();
}


inline const Foam::vector& Foam::objectMotionSolver::v0() const
{
    return motionState0_.v();
}


inline const Foam::vector& Foam::objectMotionSolver::a0() const
{
    return motionState0_.a();
}


inline const Foam::vector& Foam::objectMotionSolver::pi0() const
{
    return motionState0_.pi();
}


inline const Foam::vector& Foam::objectMotionSolver::tau0() const
{
    return motionState0_.tau();
}


inline Foam::tensor Foam::objectMotionSolver::tConstraints() const
{
    return body_.tConstraints();
}


inline Foam::tensor Foam::objectMotionSolver::rConstraints() const
{
    return body_.rConstraints();
}


inline Foam::Tuple2<Foam::tensor, Foam::vector> Foam::objectMotionSolver::rotate
(
    const tensor& Q0,
    const vector& pi,
    const scalar deltaT
) const
{
    return body_.rotate(Q0, pi, deltaT);
}


inline void Foam::objectMotionSolver::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    body_.updateAcceleration(fGlobal, tauGlobal);
}


// ************************************************************************* //
