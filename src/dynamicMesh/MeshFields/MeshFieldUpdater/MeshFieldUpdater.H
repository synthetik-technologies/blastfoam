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
#include "mapPolyMesh.H"
#include "pointMesh.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MeshField Declaration
\*---------------------------------------------------------------------------*/

template<class GeoMesh>
class MeshFieldUpdater
:
    public MeshObject
    <
        polyMesh,
        DistributeableMeshObject,
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

        virtual void updateMesh(const mapPolyMesh& map)
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

        virtual void distribute(const mapDistributePolyMesh& map);
};

#define updateFields(Type, GeoMesh)                                            \
{                                                                              \
    wordList fieldNames = mesh_.names<MeshField<Type, GeoMesh>>();             \
    forAll(fieldNames, i)                                                      \
    {                                                                          \
        MeshField<Type, GeoMesh>& f =                                          \
            mesh_.lookupObjectRef<MeshField<Type, GeoMesh>>(fieldNames[i]);    \
        mapper(f, f);                                                          \
    }                                                                          \
}

#define distributeFields(Type, GeoMesh, Distribute)                            \
{                                                                              \
    wordList fieldNames = mesh_.names<MeshField<Type, GeoMesh>>();             \
    forAll(fieldNames, i)                                                      \
    {                                                                          \
        MeshField<Type, GeoMesh>& f =                                          \
            mesh_.lookupObjectRef<MeshField<Type, GeoMesh>>(fieldNames[i]);    \
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
