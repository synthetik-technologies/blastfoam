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


Class
    Foam::MeshFieldUpdater

SourceFiles
    MeshFieldUpdater.C
    MeshFieldUpdaters.C

Description
    Parent class for updating MeshFields

\*---------------------------------------------------------------------------*/

#ifndef MeshFieldUpdater_H
#define MeshFieldUpdater_H

#include "RefineBalanceMeshObject.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MeshField Declaration
\*---------------------------------------------------------------------------*/

template<class GeoMesh>
class MeshFieldUpdater
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        RepatchableMeshObject,
        MeshFieldUpdater<GeoMesh>
    >
{
public:

    // Public typedefs
    typedef typename GeoMesh::Mesh Mesh;

    TypeName("MeshFieldUpdater");

    // Constructors

        //- Construct from IOobject
        MeshFieldUpdater(const polyMesh&);


    //- Destructor
    virtual ~MeshFieldUpdater();


    // Member Functions

        virtual bool movePoints()
        {
            return false;
        }

        virtual void topoChange(const polyTopoChangeMap& map)
        {
            NotImplemented;
        }
        virtual void mapMesh(const polyMeshMap& map)
        {
            NotImplemented;
        }

        virtual void addPatch(const label patchi)
        {}

        virtual void reorderPatches
        (
            const labelUList& newToOld,
            const bool validBoundary
        )
        {}

        virtual void distribute(const polyDistributionMap& map);
};

#define updateFields(Type, GeoMesh)                                            \
{                                                                              \
    typedef MeshField<Type, GeoMesh> FieldType;                                \
    HashTable<FieldType*> fields =                                             \
        const_cast<objectRegistry&>(mesh().thisDb()).lookupClass<FieldType>(); \
    forAllConstIter(HashTable<FieldType*>, fields, iter)                       \
    {                                                                          \
        MeshField<Type, GeoMesh>& f = *(iter());                               \
        mapper(f, f);                                                          \
    }                                                                          \
}

#define distributeFields(Type, GeoMesh, Distribute)                            \
{                                                                              \
    typedef MeshField<Type, GeoMesh> FieldType;                                \
    HashTable<FieldType*> fields =                                             \
        const_cast<objectRegistry&>(mesh().thisDb()).lookupClass<FieldType>(); \
    forAllConstIter(HashTable<FieldType*>, fields, iter)                       \
    {                                                                          \
        MeshField<Type, GeoMesh>& f = *(iter());                               \
        mdpm.Distribute(f);                                                    \
    }                                                                          \
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshFieldUpdater.C"
#endif

#endif

// ************************************************************************* //
