/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "MeshField.H"
#include "MeshFieldUpdater.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::MeshField<Type, GeoMesh>::MeshField
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    IOField<Type>(io, GeoMesh::size(mesh)),
    mesh_(mesh)
{
    MeshFieldUpdater<GeoMesh>::New(mesh_);
}


template<class Type, class GeoMesh>
Foam::MeshField<Type, GeoMesh>::MeshField
(
    const IOobject& io,
    const polyMesh& mesh,
    const Type& val
)
:
    IOField<Type>(io, GeoMesh::size(mesh)),
    mesh_(mesh)
{
    *this = val;
    MeshFieldUpdater<GeoMesh>::New(mesh_);
}


template<class Type, class GeoMesh>
Foam::MeshField<Type, GeoMesh>::MeshField
(
    const IOobject& io,
    const polyMesh& mesh,
    const Field<Type>& f
)
:
    IOField<Type>(io, f),
    mesh_(mesh)
{
    checkSize(*this);
    MeshFieldUpdater<GeoMesh>::New(mesh_);
}


template<class Type, class GeoMesh>
Foam::MeshField<Type, GeoMesh>::MeshField
(
    const IOobject& io,
    const polyMesh& mesh,
    Field<Type>&& f
)
:
    IOField<Type>(io, f),
    mesh_(mesh)
{
    checkSize(*this);
    MeshFieldUpdater<GeoMesh>::New(mesh_);
}


template<class Type, class GeoMesh>
Foam::MeshField<Type, GeoMesh>::MeshField
(
    const IOobject& io,
    const polyMesh& mesh,
    const tmp<Field<Type>>& f
)
:
    IOField<Type>(io, f),
    mesh_(mesh)
{
    checkSize(*this);
    MeshFieldUpdater<GeoMesh>::New(mesh_);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::MeshField<Type, GeoMesh>::~MeshField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
