/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
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

#include "StationaryImmersedBoundaryObject.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::
StationaryImmersedBoundaryObject
(
    const polyMesh& mesh,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(mesh, dict, stateDict),
    orientation_
    (
        dict.lookupOrDefault<tensor>("orientation", tensor::I)
      & this->shape().orientation()
    )
{}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::initialize()
{
    ImmersedType::initialize();
    this->shape_->movePoints();
    this->shape_->writeVTK();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::
~StationaryImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::Tuple2<Foam::tensor, Foam::vector>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::rotate
(
    const tensor& Q0,
    const vector& pi0,
    const scalar deltaT
) const
{
    return Tuple2<tensor, vector>(Q0, pi0);
}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    return;
}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::applyRestraints()
{
    return;
}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::status() const
{
    ImmersedType::status();
}


template<class ImmersedType>
Foam::point
Foam::StationaryImmersedBoundaryObject<ImmersedType>::centreOfRotation() const
{
    return this->initialCentreOfMass();
}


template<class ImmersedType>
Foam::point
Foam::StationaryImmersedBoundaryObject<ImmersedType>::centreOfMass() const
{
    return this->initialCentreOfMass();
}


template<class ImmersedType>
Foam::tensor
Foam::StationaryImmersedBoundaryObject<ImmersedType>::orientation() const
{
    return orientation_;
}


template<class ImmersedType>
Foam::vector
Foam::StationaryImmersedBoundaryObject<ImmersedType>::momentArm() const
{
    return Zero;
}


template<class ImmersedType>
Foam::vector
Foam::StationaryImmersedBoundaryObject<ImmersedType>::v() const
{
    return Zero;
}


template<class ImmersedType>
Foam::vector
Foam::StationaryImmersedBoundaryObject<ImmersedType>::velocity
(
    const point& pt
) const
{
    return Zero;
}


template<class ImmersedType>
Foam::tmp<Foam::vectorField>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::velocity
(
    const pointField& points
) const
{
    return tmp<vectorField>(new vectorField(points.size(), Zero));
}


template<class ImmersedType>
Foam::vector
Foam::StationaryImmersedBoundaryObject<ImmersedType>::omega() const
{
    return Zero;
}


template<class ImmersedType>
Foam::diagTensor
Foam::StationaryImmersedBoundaryObject<ImmersedType>::momentOfInertia() const
{
    return diagTensor(great, great, great);
}


template<class ImmersedType>
Foam::point
Foam::StationaryImmersedBoundaryObject<ImmersedType>::transform
(
    const point& initialPoint
) const
{
    return centreOfRotation() + (orientation() & initialPoint);
}


template<class ImmersedType>
Foam::tmp<Foam::pointField>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::transform
(
    const pointField& initialPoints
) const
{
    return centreOfRotation() + (orientation() & initialPoints);
}


template<class ImmersedType>
Foam::point Foam::StationaryImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const point& pt
) const
{
    return (orientation().inv() & (pt - centreOfRotation()));
}


template<class ImmersedType>
Foam::tmp<Foam::pointField> Foam::StationaryImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const pointField& points
) const
{
    return (orientation().inv() & (points - centreOfRotation()));
}


template<class ImmersedType>
Foam::tensor
Foam::StationaryImmersedBoundaryObject<ImmersedType>::rConstraints() const
{
    return tensor::I;
}


template<class ImmersedType>
Foam::tensor
Foam::StationaryImmersedBoundaryObject<ImmersedType>::tConstraints() const
{
    return tensor::I;
}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::movePoints()
{
    ImmersedType::movePoints();
}


template<class ImmersedType>
bool Foam::StationaryImmersedBoundaryObject<ImmersedType>::read
(
    const dictionary& dict
)
{
    return ImmersedType::read(dict);
}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
