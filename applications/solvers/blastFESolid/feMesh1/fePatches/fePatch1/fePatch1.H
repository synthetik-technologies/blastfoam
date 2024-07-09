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
    Foam::fePatch1

Description
    Basic fePatch1 represents a set of points from the mesh.

SourceFiles
    fePatch1.C

\*---------------------------------------------------------------------------*/

#ifndef fePatch1_H
#define fePatch1_H

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

class feBoundaryMesh1;
class pointConstraint;
class PstreamBuffers;

/*---------------------------------------------------------------------------*\
                         Class fePatch1 Declaration
\*---------------------------------------------------------------------------*/

class fePatch1
{
    // Private Data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;

        //- Reference to boundary mesh
        const feBoundaryMesh1& boundaryMesh_;


// Private Data

    // Demand driven private data

        //- Elements
        mutable List<element>* elementsPtr_;

        //- Local element indicies
        mutable List<labelList>* localElementsPtr_;

        //- Local nodes
        mutable pointField* localNodesPtr_;

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

        //- Calculate local elements
        void calcLocalElements() const;

        //- Calculate local nodes
        void calcLocalNodes() const;

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

        // The fePatch1 geometry initialisation is called by pointBoundaryMesh
        friend class feBoundaryMesh1;

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&)
        {
            deleteDemandDrivenData(elementsPtr_);
            deleteDemandDrivenData(localElementsPtr_);
            deleteDemandDrivenData(localNodesPtr_);
            deleteDemandDrivenData(shapesPtr_);
            deleteDemandDrivenData(dshapesPtr_);
            deleteDemandDrivenData(BsPtr_);
            deleteDemandDrivenData(JsPtr_);
            deleteDemandDrivenData(invJsPtr_);
            deleteDemandDrivenData(WsPtr_);
            deleteDemandDrivenData(nodeNormalsPtr_);
        }

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {
            deleteDemandDrivenData(BsPtr_);
            deleteDemandDrivenData(JsPtr_);
            deleteDemandDrivenData(invJsPtr_);
            deleteDemandDrivenData(WsPtr_);
            deleteDemandDrivenData(nodeNormalsPtr_);
        }

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&)
        {
            deleteDemandDrivenData(elementsPtr_);
            deleteDemandDrivenData(shapesPtr_);
            deleteDemandDrivenData(dshapesPtr_);
            deleteDemandDrivenData(BsPtr_);
            deleteDemandDrivenData(JsPtr_);
            deleteDemandDrivenData(invJsPtr_);
            deleteDemandDrivenData(WsPtr_);
            deleteDemandDrivenData(nodeNormalsPtr_);
        }

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&)
        {}


public:

    //- Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fePatch1,
            polyPatch,
            (const polyPatch& patch, const feBoundaryMesh1& bm),
            (patch, bm)
        );

    // Constructor

        fePatch1(const polyPatch& pp, const feBoundaryMesh1& bm);

        //- Disallow default bitwise copy construction
        fePatch1(const fePatch1&) = delete;


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<fePatch1> New
        (
            const polyPatch&,
            const feBoundaryMesh1&
        );


    //- Destructor
    virtual ~fePatch1();


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
            const feBoundaryMesh1& boundaryMesh() const
            {
                return boundaryMesh_;
            }


        // Access functions for Topological data

            //- Return number of points supporting patch faces
            label nNodes() const;

            //- Access elements
            const List<element>& elements() const;

            //- Return local elements
            const List<labelList>& localElements() const;

            //- Return local elements
            const pointField& localNodes() const;


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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fePatch1&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
