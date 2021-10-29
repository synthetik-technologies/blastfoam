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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ImmersedType>
inline Foam::tensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::rotationTensorX
(
    scalar phi
) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


template<class ImmersedType>
inline Foam::tensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::rotationTensorY
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


template<class ImmersedType>
inline Foam::tensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::rotationTensorZ
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


template<class ImmersedType>
inline const Foam::PtrList<Foam::immersedBoundaryObjectRestraint>&
Foam::MovingImmersedBoundaryObject<ImmersedType>::restraints() const
{
    return restraints_;
}


template<class ImmersedType>
inline const Foam::PtrList<Foam::immersedBoundaryObjectConstraint>&
Foam::MovingImmersedBoundaryObject<ImmersedType>::constraints() const
{
    return constraints_;
}


template<class ImmersedType>
inline const Foam::tensor&
Foam::MovingImmersedBoundaryObject<ImmersedType>::Q() const
{
    return motionState_.Q();
}


template<class ImmersedType>
inline const Foam::vector&
Foam::MovingImmersedBoundaryObject<ImmersedType>::a() const
{
    return motionState_.a();
}


template<class ImmersedType>
inline const Foam::vector&
Foam::MovingImmersedBoundaryObject<ImmersedType>::pi() const
{
    return motionState_.pi();
}


template<class ImmersedType>
inline const Foam::vector&
Foam::MovingImmersedBoundaryObject<ImmersedType>::tau() const
{
    return motionState_.tau();
}


template<class ImmersedType>
inline Foam::tensor& Foam::MovingImmersedBoundaryObject<ImmersedType>::Q()
{
    return motionState_.Q();
}


template<class ImmersedType>
inline Foam::vector& Foam::MovingImmersedBoundaryObject<ImmersedType>::v()
{
    return motionState_.v();
}


template<class ImmersedType>
inline Foam::vector& Foam::MovingImmersedBoundaryObject<ImmersedType>::a()
{
    return motionState_.a();
}


template<class ImmersedType>
inline Foam::vector& Foam::MovingImmersedBoundaryObject<ImmersedType>::pi()
{
    return motionState_.pi();
}


template<class ImmersedType>
inline Foam::vector& Foam::MovingImmersedBoundaryObject<ImmersedType>::tau()
{
    return motionState_.tau();
}

template<class ImmersedType>
inline Foam::vector&
Foam::MovingImmersedBoundaryObject<ImmersedType>::centreOfRotation()
{
    return motionState_.centreOfRotation();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
inline const Foam::movingImmersedBoundaryObjectState&
Foam::MovingImmersedBoundaryObject<ImmersedType>::state() const
{
    return motionState_;
}


template<class ImmersedType>
inline void Foam::MovingImmersedBoundaryObject<ImmersedType>::newTime()
{
    motionState0_ = motionState_;
}


template<class ImmersedType>
inline const Foam::pointField&
Foam::MovingImmersedBoundaryObject<ImmersedType>::points0() const
{
    return points0_;
}

// ************************************************************************* //
