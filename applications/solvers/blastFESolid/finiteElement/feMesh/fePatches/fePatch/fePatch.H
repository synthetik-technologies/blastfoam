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
    Foam::fePatch

Description
    Basic fePatch represents a set of points from the mesh.

SourceFiles
    fePatch.C

\*---------------------------------------------------------------------------*/

#ifndef fePatch_H
#define fePatch_H

#include "polyPatch.H"
#include "Time.H"
#include "element.H"
#include "labelList.H"
#include "vectorField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class feBoundaryMesh;
class pointConstraint;
class PstreamBuffers;

/*---------------------------------------------------------------------------*\
                         Class fePatch Declaration
\*---------------------------------------------------------------------------*/

class fePatch
{
    // Private Data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;

        //- Reference to boundary mesh
        const feBoundaryMesh& boundaryMesh_;


// Private Data

    // Demand driven private data

        //- Elements
        mutable List<element>* elementsPtr_;

        //- Boundary point labels, addressing into local point list
        mutable labelList* boundaryNodesPtr_;

        //- Labels of mesh nodes
        mutable labelList* meshNodesPtr_;

        //- Saved shape functions
        mutable List<List<scalarList>>* shapesPtr_;

        //- Saved derivative of shape functions
        mutable List<List<scalarRectangularMatrix>>* dshapesPtr_;

        //- Saved physcial gradient of shape functions
        mutable List<List<scalarRectangularMatrix>>* BsPtr_;

        //- Saved Jacobains
        mutable List<List<tensor>>* JsPtr_;

        //- Saved inverse Jacobains
        mutable List<List<tensor>>* invJsPtr_;

        //- Saved determinates of Jacobians
        mutable List<List<scalar>>* WsPtr_;

        //- Node unit normals
        mutable List<List<vector>>* nodeNormalsPtr_;


    // Private Member Functions

        //- Calculate elements
        void calcElements() const;

        //- Calculate mesh addressing
        void calcMeshData() const;

        //- Calculate shape functions
        void calcShapes() const;

        //- Calculate derivative of shape functions
        void calcDShapes() const;

        //- Calculate Jacobians
        void calcJs() const;

        //- Calculate inverse Jacobians
        void calcInvJs() const;

        //- Calculate physical gradients of shape functions
        void calcBs() const;

        //- Calculate weights
        void calcWs() const;

        //- Calculate unit node normals
        void calcNodeNormals() const;


protected:

    // Protected Member Functions

        // The fePatch geometry initialisation is called by pointBoundaryMesh
        friend class feBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&)
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&)
        {}

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&)
        {}

        //- Clear addressing
        void clearAddressing() const;

        //- Clear geometric data
        void clearGeometry() const;


public:

    //- Runtime type information
    TypeName("fePatch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fePatch,
            polyPatch,
            (const polyPatch& patch, const feBoundaryMesh& bm),
            (patch, bm)
        );

    // Constructor

        fePatch(const polyPatch& pp, const feBoundaryMesh& bm);

        //- Disallow default bitwise copy construction
        fePatch(const fePatch&) = delete;


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<fePatch> New
        (
            const polyPatch&,
            const feBoundaryMesh&
        );


    //- Destructor
    virtual ~fePatch();


    // Member Functions

        // Access

            //- Return the polyPatch
            const polyPatch& patch() const
            {
                return polyPatch_;
            }

            //- Return name
            virtual const word& name() const
            {
                return polyPatch_.name();
            }

            //- Return size
            virtual label size() const
            {
                return polyPatch_.size();
            }

            //- Return true if this patch is coupled
            virtual bool coupled() const
            {
                return polyPatch_.coupled();
            }

            //- Return the index of this patch in the fvBoundaryMesh
            label index() const
            {
                return polyPatch_.index();
            }

            //- Return boundaryMesh reference
            const feBoundaryMesh& boundaryMesh() const
            {
                return boundaryMesh_;
            }


        // Access functions for Topological data

            //- Return number of points supporting patch faces
            label nNodes() const;

            //- Access elements
            const List<element>& elements() const;


        // Addressing into mesh

            //- Return labelList of mesh nodes in patch. They are constructed
            //  walking through the faces in incremental order and not sorted
            //  anymore.
            const labelList& meshNodes() const;

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

            //- Accumulate the effect of constraint direction of this patch
            virtual void applyConstraint
            (
                const label nodei,
                pointConstraint&
            ) const
            {}

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fePatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
