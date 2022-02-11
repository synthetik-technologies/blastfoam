/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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
    Foam::BalanceMeshObject
    Foam::RefineMeshObject
    Foam::BalanceMeshObject

Description
    Object used to lookup after refinement and balancing

SourceFiles
    RefineBalanceMeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef RefineBalanceMeshObject_H
#define RefineBalanceMeshObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh.H"
#include "MeshObject.H"
#include "mapDistributePolyMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class DistributeableMeshObject Declaration
\*---------------------------------------------------------------------------*/

class blastMeshObject
{
public:
    template<class Mesh>
    static void distribute
    (
        objectRegistry& obr,
        const mapDistributePolyMesh& map
    );
};

template<class Mesh>
class DistributeableMeshObject
:
    public PatchMeshObject<Mesh>
{
public:

    DistributeableMeshObject(const word& typeName, const objectRegistry& obr)
    :
        PatchMeshObject<Mesh>(typeName, obr)
    {}

    DistributeableMeshObject(const word& typeName, const IOobject& io)
    :
        PatchMeshObject<Mesh>(typeName, io)
    {}

    virtual void distribute(const mapDistributePolyMesh& map) = 0;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "RefineBalanceMeshObject.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
