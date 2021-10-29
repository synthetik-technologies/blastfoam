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

#include "DeformableImmersedMesh.H"
#include "dynamicBlastFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::DeformableImmersedMesh<ImmersedType>::
DeformableImmersedMesh
(
    const polyMesh& mesh,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(mesh, dict, stateDict),
    fvMesh_
    (
        immersedShape::castShapeType<immersedFvMesh>
        (
            this->shape()
        ).immersedDyMesh()
    ),
    solidPtr_(),
    mapperPtr_(nullptr),
    patchName_(dict.lookupOrDefault("patchName", word::null))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::DeformableImmersedMesh<ImmersedType>::~DeformableImmersedMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void Foam::DeformableImmersedMesh<ImmersedType>::solve()
{
    solidPtr_->evolve();
}

template<class ImmersedType>
const Foam::fvMesh&
Foam::DeformableImmersedMesh<ImmersedType>::immersedMesh() const
{
    return fvMesh_;
}


template<class ImmersedType>
const Foam::immersedMeshMapper*
Foam::DeformableImmersedMesh<ImmersedType>::mapper() const
{
    return mapperPtr_;
}


template<class ImmersedType>
void Foam::DeformableImmersedMesh<ImmersedType>::movePoints()
{
    if (!solidPtr_->movingMesh())
    {
        solidPtr_->moveMesh(points0Ptr_(), solidPtr_->D(), solidPtr_->pointD());
    }

    this->shape_->movePoints();

    this->clearOut();

    if (!solidPtr_->movingMesh())
    {
        fvMesh_.movePoints(points0Ptr_());
        fvMesh_.V00();
        fvMesh_.moving(false);
        fvMesh_.setPhi().writeOpt() = IOobject::NO_WRITE;
    }


}


template<class ImmersedType>
void Foam::DeformableImmersedMesh<ImmersedType>::initialize()
{
    mapperPtr_ = new immersedMeshMapper(fvMesh_, *this, patchName_);
    patchName_ = fvMesh_.boundaryMesh()[mapperPtr_->interfaceIndex()].name();
    solidPtr_ = solidModel::New(fvMesh_);

    if (!solidPtr_->movingMesh())
    {
        points0Ptr_.set(new pointField(fvMesh_.points()));
        mapperPtr_->setDisplacement
        (
            solidPtr_->pointD(),
            solidPtr_->D()
        );
    }
}


template<class ImmersedType>
Foam::scalar Foam::DeformableImmersedMesh<ImmersedType>::CoNum() const
{
   return solidPtr_->CoNum();
}


template<class ImmersedType>
Foam::scalar Foam::DeformableImmersedMesh<ImmersedType>::maxCoNum() const
{
    return solidPtr_->maxCoNum();
}


template<class ImmersedType>
void Foam::DeformableImmersedMesh<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::DeformableImmersedMesh<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
