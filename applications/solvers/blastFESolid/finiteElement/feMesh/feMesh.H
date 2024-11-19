/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::feMesh

Description
    Mesh data needed to do the Finite Element discretisation.

    NOTE ON USAGE:
    feMesh contains all the topological and geometric information
    related to the mesh.  It is also responsible for keeping the data
    up-to-date.  This is done by deleting the cell volume, face area,
    cell/face centre, addressing and other derived information as
    required and recalculating it as necessary.  The feMesh therefore
    reserves the right to delete the derived information upon every
    topological (mesh refinement/morphing) or geometric change (mesh
    motion).  It is therefore unsafe to keep local references to the
    derived data outside of the time loop.

SourceFiles
    feMesh.C
    feMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef feMesh_H
#define feMesh_H

#include "GeoMesh.H"
#include "MeshObject.H"
#include "polyMesh.H"
#include "feBoundaryMesh.H"
#include "element.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class feMesh Declaration
\*---------------------------------------------------------------------------*/

class feMesh
:
    public objectRegistry,
    public MeshObject<polyMesh, PatchMeshObject, feMesh>,
    public GeoMesh<polyMesh>
{
    // Private Data

        //- Mesh order
        label order_;

        //- Boundary mesh
        feBoundaryMesh boundary_;

    // Demand-driven data

        //- Volume elements
        mutable List<element>* elementsPtr_;

        //- Nodes, only set if order > 1
        mutable pointIOField* nodesPtr_;

        //- Edge-node indexing
        mutable labelListList* edgeNodesPtr_;

        //- Face-node indexing
        mutable labelListList* faceNodesPtr_;

        //- Saved shape functions
        //  elements x integration points
        mutable List<List<scalarList>>* shapesPtr_;

        //- Saved derivative of shape functions
        //  elements x integration points
        mutable List<List<scalarRectangularMatrix>>* dshapesPtr_;

        //- Saved physcial gradient of shape functions
        //  elements x integration points
        mutable List<List<scalarRectangularMatrix>>* BsPtr_;

        //- Saved Jacobains
        //  elements x integration points
        mutable List<List<tensor>>* JsPtr_;

        //- Saved inverse Jacobains
        //  elements x integration points
        mutable List<List<tensor>>* invJsPtr_;

        //- Saved determinates of Jacobians
        //  elements x integration points
        mutable List<List<scalar>>* WsPtr_;

        //- Number of integration points
        mutable label nIp_;

        //- Element integration points
        mutable List<labelList> ipLabelsPtr_;


    // Private Member Functions

        // Storage management

            //- Clear geometry but not the old-time cell volumes
            void clearGeom();

            //- Clear addressing
            void clearAddressing(const bool isMeshUpdate = false);

            //- Clear Geometry and addressing
            void clearOut();


       // Make geometric data
            void calcElements() const;
            void calcIpLabels() const;
            void calcShapes() const;
            void calcDShapes() const;
            void calcBs() const;

            void calcJs() const;
            void calcInvJs() const;
            void calcWs() const;

        //- Map fields
        void mapFields(const mapPolyMesh& mpm);

public:

    // Public Typedefs

        typedef feMesh Mesh;
        typedef feBoundaryMesh BoundaryMesh;


    // Declare name of the class and its debug switch
    ClassName("feMesh");


    // Constructors

        //- Construct from IOobject
        explicit feMesh(const label order, const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        feMesh(const feMesh&);


    //- Destructor
    virtual ~feMesh();


    // Member Functions

        //- Return the time
        const Time& time() const
        {
            return objectRegistry::time();
        }

        //- Return number of points
        label size() const
        {
            return size(*this);
        }

        //- Return number of points
        static label size(const Mesh& mesh)
        {
            return mesh.GeoMesh<polyMesh>::mesh_.nCells();
        }

        //- Return reference to boundary mesh
        const feBoundaryMesh& boundary() const
        {
            return boundary_;
        }

        //- Return parallel info
        const globalMeshData& globalData() const
        {
            return GeoMesh<polyMesh>::mesh_.globalData();
        }

        //- Return database. For now is its polyMesh.
        const objectRegistry& thisDb() const
        {
            return GeoMesh<polyMesh>::mesh_.thisDb();
        }

        //- Reset for changed polyMesh
        void reset(const bool validBoundary);


        // Access functions for Topological data

            //- Return mesh order
            label order() const
            {
                return order_;
            }

            //- Return number of points supporting patch faces
            label nNodes() const;

            //- Number of integration points
            label nIp() const

            //- Access elements
            const List<element>& elements() const;

            //- Access element integration point labels
            const List<labelList>& ipLabels() const;

            //- Return nodes
            const pointField& nodes() const;

            //- Access edge-node addressing
            const List<labelList>& edgeNodes() const;

            //- Access face-node addressing
            const List<labelList>& faceNodes() const;

            //- Access face-node addressing
            const List<labelList>& faceNodes() const;


        // Addressing into mesh

            //- Return shape functions
            const List<List<scalarList>>& shapes() const;

            //- Return derivative of shape functions
            const List<List<scalarRectangularMatrix>>& dshapes() const;

            //- Return physical gradient of shape functions
            const List<List<scalarRectangularMatrix>>& Bs() const;

            //- Return Jacobians
            const List<List<tensor>>& Js() const;

            //- Return inverse Jacobians
            const List<List<tensor>>& invJs() const;

            //- Return weights
            const List<List<scalar>>& Ws() const;

            //- Return node normals for patch
            const List<List<vector>>& nodeNormals() const;


        // Mesh callbacks

            //- Move points
            virtual bool movePoints();

            //- Update the mesh corresponding to given map
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Reordered/removed trailing patches. If validBoundary call is
            //  parallel synced and all add the same patch with same settings
            virtual void reorderPatches
            (
                const labelUList& newToOld,
                const bool validBoundary
            );

            //- Inserted patch at patchi
            virtual void addPatch(const label patchi);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const feMesh&);

        bool operator!=(const feMesh&) const;
        bool operator==(const feMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
