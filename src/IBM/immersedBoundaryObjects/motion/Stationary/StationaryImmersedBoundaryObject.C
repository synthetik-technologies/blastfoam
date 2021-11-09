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
    const polyPatch& patch,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(patch, dict, stateDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::
~StationaryImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void Foam::StationaryImmersedBoundaryObject<ImmersedType>::status
(
    const bool print
) const
{
    ImmersedType::status(false);
}


template<class ImmersedType>
Foam::point
Foam::StationaryImmersedBoundaryObject<ImmersedType>::centre() const
{
    return this->shape_->centre();
}


template<class ImmersedType>
Foam::scalar
Foam::StationaryImmersedBoundaryObject<ImmersedType>::mass() const
{
    return great;
}


template<class ImmersedType>
Foam::diagTensor
Foam::StationaryImmersedBoundaryObject<ImmersedType>::momentOfInertia() const
{
    return diagTensor(great, great, great);
}


template<class ImmersedType>
Foam::vector
Foam::StationaryImmersedBoundaryObject<ImmersedType>::v
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
Foam::tmp<Foam::vectorField>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::velocity() const
{
    return tmp<vectorField>(new vectorField(this->size(), Zero));
}


template<class ImmersedType>
Foam::tensor
Foam::StationaryImmersedBoundaryObject<ImmersedType>::orientation() const
{
    return this->shape_->orientation();
}


template<class ImmersedType>
Foam::point
Foam::StationaryImmersedBoundaryObject<ImmersedType>::transform
(
    const point& initialPoint
) const
{
    return centre() + (orientation() & initialPoint);
}


template<class ImmersedType>
Foam::tmp<Foam::pointField>
Foam::StationaryImmersedBoundaryObject<ImmersedType>::transform
(
    const pointField& initialPoints
) const
{
    return centre() + (orientation() & initialPoints);
}


template<class ImmersedType>
Foam::point Foam::StationaryImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const point& pt
) const
{
    return (orientation().T() & (pt - centre()));
}


template<class ImmersedType>
Foam::tmp<Foam::pointField> Foam::StationaryImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const pointField& points
) const
{
    return (orientation().T() & (points - centre()));
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
