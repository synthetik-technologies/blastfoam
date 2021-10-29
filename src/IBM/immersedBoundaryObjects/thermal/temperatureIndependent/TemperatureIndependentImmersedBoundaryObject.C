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

#include "TemperatureIndependentImmersedBoundaryObject.H"
#include "fvm.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::
TemperatureIndependentImmersedBoundaryObject
(
    const polyMesh& mesh,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(mesh, dict, stateDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::
~TemperatureIndependentImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::status() const
{
    ImmersedType::status();
}


template<class ImmersedType>
void Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::movePoints()
{
    ImmersedType::movePoints();
}


template<class ImmersedType>
bool Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::read
(
    const dictionary& dict
)
{
    return ImmersedType::read(dict);
}


template<class ImmersedType>
void Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::TemperatureIndependentImmersedBoundaryObject<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
