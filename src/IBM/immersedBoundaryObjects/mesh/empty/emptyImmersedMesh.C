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

#include "NondeformableImmersedMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::NondeformableImmersedMesh<ImmersedType>::
NondeformableImmersedMesh
(
    const polyMesh& mesh,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(mesh, dict, stateDict),
    fvMesh_
    (
        IOobject
        (
            this->name(),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    ),
    mapperPtr_(nullptr),
    points0_(fvMesh_.points())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::NondeformableImmersedMesh<ImmersedType>::
~NondeformableImmersedMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
const Foam::fvMesh*
Foam::NondeformableImmersedMesh<ImmersedType>::immersedFvMesh() const
{
    return &fvMesh_;
}


template<class ImmersedType>
const Foam::immersedMeshMapper*
Foam::NondeformableImmersedMesh<ImmersedType>::mapper() const
{
    return mapperPtr_;
}


template<class ImmersedType>
void Foam::NondeformableImmersedMesh<ImmersedType>::movePoints()
{
    fvMesh_.movePoints(this->transform(points0_));
    fvMesh_.moving(false);
    ImmersedType::movePoints();
}


template<class ImmersedType>
void Foam::NondeformableImmersedMesh<ImmersedType>::initialize()
{
    mapperPtr_ = new immersedMeshMapper(fvMesh_, *this);
}


template<class ImmersedType>
void Foam::NondeformableImmersedMesh<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::NondeformableImmersedMesh<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
