/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License

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
    Foam::globalPolyBoundaryMesh

Description
    Mesh object used to collect all mapped patches so that they can be cleared


Note

SourceFiles
    globalPolyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef globalPolyBoundaryMesh_H
#define globalPolyBoundaryMesh_H

#include "globalPolyPatch.H"
#include "regIOobject.H"
#include "HashPtrTable.H"
#include "MeshObject.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointPatch;
class coupledGlobalPolyPatch;

class globalPolyBoundaryMesh;
typedef MeshObject
<
    polyMesh,
    PatchMeshObject,
    globalPolyBoundaryMesh
> GlobalPolyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                       Class globalPolyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class globalPolyBoundaryMesh
:
    public GlobalPolyBoundaryMesh
{
// Private member data

    const HashTable<dictionary> interfacesDicts_;

    // Table of mapped patch selectors
    mutable HashPtrTable<globalPolyPatch> patches_;


public:

    //- Runtime type information
    TypeName("globalPolyBoundaryMesh");


    // Constructors

        //- Construct from patch
        globalPolyBoundaryMesh(const polyMesh&);


    //- Destructor
    virtual ~globalPolyBoundaryMesh();


    // Move points
    virtual bool movePoints();

    //- Update mesh
    virtual void updateMesh(const mapPolyMesh& mpm);

    //- Reorder patches
    virtual void reorderPatches
    (
        const labelUList& newToOld,
        const bool validBoundary
    );

    //- Inserted patch at patchi
    virtual void addPatch(const label patchi);

    // Operators

        //- Return the globalPolyPatch for the given polyPatch
        const globalPolyPatch& operator[](const word&) const;

        //- Return the globalPolyPatch for the given polyPatch
        const globalPolyPatch& operator[](const polyPatch&) const;

        //- Return the globalPolyPatch for the given pointPatch
        const globalPolyPatch& operator[](const pointPatch&) const;

        //- Return the coupledGlobalPolyPatch for the given polyPatch
        const coupledGlobalPolyPatch& operator()(const word&) const;

        //- Return the coupledGlobalPolyPatch for the given polyPatch
        const coupledGlobalPolyPatch& operator()(const polyPatch&) const;

        //- Return the coupledGlobalPolyPatch for the given pointPatch
        const coupledGlobalPolyPatch& operator()(const pointPatch&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
