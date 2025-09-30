/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::feMesh1

Description
    Mesh data needed to do the Finite Element discretisation. Only 1st
    order geometry is used


SourceFiles
    feMesh1.C
    feMesh1Geometry.C

\*---------------------------------------------------------------------------*/

#ifndef feMesh1_H
#define feMesh1_H

#include "pointMesh.H"
#include "MeshObjects.H"
#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "feBoundaryMesh1.H"
#include "element.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class feMesh1 Declaration
\*---------------------------------------------------------------------------*/

class feMesh1
:
    public DemandDrivenMeshObject<polyMesh, RepatchableMeshObject, feMesh1>
{
    // Private Data

        //- Reference to pointMesh
        const pointMesh& pointMesh_;

        //- Integration order
        const label intOrder_;

        //- Shell integration order
        const label shellIntOrder_;

        //- Is this a shell mesh
        const bool shell_;

        //- Boundary mesh
        feBoundaryMesh1 boundary_;

    // Demand-driven data

        //- Volume elements
        mutable List<element>* elementsPtr_;

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

        //- Weight field (row summation)
        mutable pointScalarField* WPtr_;

        //- Number of integration points
        mutable label nIp_;

        //- Element integration points
        mutable List<labelList>* ipLabelsPtr_;


    // Private Member Functions

        // Storage management

            //- Clear reference geometry
            void clearRefGeom();

            //- Clear physical geometry
            void clearPhysGeom();

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
            void calcW() const;

        //- Map fields
        void mapFields(const polyTopoChangeMap& map);

public:

    // Public Typedefs

        typedef feMesh1 Mesh;
        typedef feBoundaryMesh1 BoundaryMesh;


    // Declare name of the class and its debug switch
    ClassName("feMesh1");


    // Constructors

        //- Construct from IOobject
        explicit feMesh1
        (
            const polyMesh& mesh,
            const label intOrder
        );

        explicit feMesh1
        (
            const polyMesh& mesh,
            const label orderRS,
            const label orderT
        );

        //- Disallow default bitwise copy construction
        feMesh1(const feMesh1&);


    //- Destructor
    virtual ~feMesh1();


    // Member Functions

        //- Return reference to point mesh
        const pointMesh& pMesh() const
        {
            return pointMesh_;
        }

        //- Return reference to point boundary mesh
        const pointBoundaryMesh& pointBoundary() const
        {
            return pointMesh_.boundary();
        }

        //- Return reference to boundary mesh
        const feBoundaryMesh1& boundary() const
        {
            return boundary_;
        }

        //- Reset for changed polyMesh
        void reset(const bool validBoundary);


        // Access functions for Topological data

            //- Return integration order
            label intOrder() const
            {
                return intOrder_;
            }

            //- Return number of points supporting patch faces
            label nNodes() const;

            //- Number of integration points
            label nIp() const;

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

            //- Return weight field (row summation)
            const pointScalarField& W() const;


        // Mesh callbacks

            //- Move points
            virtual bool movePoints();

            //- Update the mesh corresponding to given map
            virtual void distribute(const polyDistributionMap& map);
            virtual void topoChange(const polyTopoChangeMap& map);
            virtual void mapMesh(const polyMeshMap& map);

            //- Reordered/removed trailing patches. If validBoundary call is
            //  parallel synced and all add the same patch with same settings
            virtual void reorderPatches
            (
                const labelUList& newToOld,
                const bool validBoundary
            );

            //- Inserted patch at patchi
            virtual void addPatch(const label patchi);


        //- Synchronize data between processor boundaries
        template<class Type>
        void pushUntransformedData(Field<Type>& pointData) const
        {
            const globalMeshData& gmd = mesh().globalData();

            // Transfer onto coupled patch

            const indirectPrimitivePatch& cpp = gmd.coupledPatch();
            const labelList& meshPoints = cpp.meshPoints();

            const distributionMap& slavesMap = gmd.globalCoPointSlavesMap();
            const labelListList& slaves = gmd.globalCoPointSlaves();

            List<Type> elems(slavesMap.constructSize());
            forAll(meshPoints, i)
            {
                elems[i] = pointData[meshPoints[i]];
            }

            // Combine master data with slave data
            forAll(slaves, i)
            {
                const labelList& slavePoints = slaves[i];

                // Copy master data to slave slots
                forAll(slavePoints, j)
                {
                    elems[slavePoints[j]] = elems[i];
                }
            }

            // Push slave-slot data back to slaves
            slavesMap.reverseDistribute(elems.size(), elems, false);

            // Extract back onto mesh
            forAll(meshPoints, i)
            {
                pointData[meshPoints[i]] = elems[i];
            }
        }



    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const feMesh1&);

        bool operator!=(const feMesh1&) const;
        bool operator==(const feMesh1&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
