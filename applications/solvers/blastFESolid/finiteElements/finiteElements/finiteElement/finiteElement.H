/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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
    Foam::finiteElement

Description
    Base class for a finite element

SourceFiles
    finiteElementI.H
    finiteElement.C

\*---------------------------------------------------------------------------*/

#ifndef finiteElement_H
#define finiteElement_H

#include "List.H"
#include "integrationRule.H"
#include "scalarMatrices.H"
#include "HashPtrTable.H"
#include "labelPair.H"
#include "runTimeSelectionTables.H"
#include "runTimeSelectionMaps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class finiteElement Declaration
\*---------------------------------------------------------------------------*/

class finiteElement
{
public:

    //- Stored elements hashed by (ElementType, order)
    typedef HashPtrTable
    <
        finiteElement,
        Tuple2<word, label>,
        Tuple2<word, label>::Hash<>
    > finiteElementTable;
    static finiteElementTable finiteElements;

protected:

    //- Order of the element
    label order_;

    //- Integration rule
    integrationRule ir_;

    //- Local nodes in reference space
    pointField nodes_;

    //- Type of the element
    ElementType::Type type_;

    //- VTK connectivity
    mutable labelList vtkConnectivity_;

    //- Uniform points
    mutable List<point> P_;

    //- Uniform edges
    mutable List<labelPair> E_;

    //- Uniform geometry
    mutable List<labelList> G_;

    //- Binomial coefficients
    static List<List<label>> Binomials;
    static const List<label>& Binomial(const label p);

    //- Binomial coefficients
    static void calcBinomialTerms
    (
        const label p,
        const scalar x,
        const scalar y,
        List<label>& u
    );

    static void calcBinomialTerms
    (
        const label p,
        const scalar x,
        const scalar y,
        List<label>& u,
        List<label>& du
    );

    //- Return the map for line/square/cube elements
    static labelList dofMap(const label nDim, const label o);


    // VTK write functions


        // From MFEM https://github.com/mfem/mfem/blob/master/mesh/vtk.hpp

        // Return the VTK node index of the barycentric point b in a
        ///triangle with refinement level.
        static label barycentricToVTKTriangle(label* b, label ref);


        static label barycentricToVTKTetra(label* b, label ref);
        static label vtkTriangleDOFOffset
        (
            const label ref,
            const label i,
            const label j
        );
        static label cartesianToVTKPrism
        (
            const label i,
            const label j,
            const label k,
            const label ref
        );
        static label cartesianToVTKTensor
        (
            const label idx_in,
            const label ref,
            const ElementType::Type geom
        );

        // Create the VTK element connectivity array for a given element
        // geometry and refinement level.
        static labelList vtkElementConnectivity
        (
            const ElementType::Type geom,
            const label ref
        );

        //- Uniform points
        static List<scalar> linspace
        (
            const label n,
            const scalar l = 0,
            const scalar u = 1
        )
        {
            scalarList pts(n);
            scalar dx = (u - l)/scalar(n - 1);
            forAll(pts, i)
            {
                pts[i] = scalar(i)*dx;
            }
            return pts;
        }

        //- Create uniform geometry
        virtual void createUniformGeometry() const = 0;

public:

    //- Runtime type information
    TypeName("finiteElement")

    // Declare runtime constructor selection tables

        //- Look up by type, given an order
        declareRunTimeSelectionTable
        (
            autoPtr,
            finiteElement,
            type,
            (
                const label order
            ),
            (order)
        );

        //- Lookup only by name
        declareRunTimeSelectionTable
        (
            autoPtr,
            finiteElement,
            typeOrder,
            (),
            ()
        );

        //- Lookup by vtk index
        declareRunTimeSelectionMap
        (
            autoPtr,
            finiteElement,
            vtk,
            (),
            ()
        );

        //- Lookup by gmsh index
        declareRunTimeSelectionMap
        (
            autoPtr,
            finiteElement,
            msh,
            (),
            ()
        );

    // Constructors

        //- Construct from order, size and element type
        finiteElement
        (
            const label order,
            const label size,
            const ElementType::Type type
        );

    //- Destructor
    virtual ~finiteElement();


    // Member functions

        // Access

            //- Return the element type
            inline ElementType::Type elemType() const;

            //- Return the element order
            inline label order() const;

            //- Return the integration rule
            inline const integrationRule& ir() const;

            // Return integration rule for a given order
            virtual const integrationRule& ir
            (
                const label orderi,
                const label orderj,
                const label orderk
            ) const = 0;

            //- Return the reference nodes
            inline const pointField& nodes() const;

            //- Return the number of dimensions
            virtual label nDims() const = 0;

            //- Return the number of vertices
            virtual label nVertices() const = 0;

            // Return number of edges
            virtual label nEdges() const = 0;

            //- Return the number of faces
            virtual label nFaces() const = 0;

            //- Return the number of nodes
            virtual label nNodes() const = 0;

            //- Return the number of edge nodes (neglecting start/end)
            virtual label nEdgeNodes() const = 0;

            //- Return the number of nodes for a face, neglecting edges
            virtual label nFaceNodes(const label fi) const = 0;

            //- Return the number of internal nodes
            virtual label nInternalNodes() const = 0;

            //- VTK index
            virtual label vtkIndex() const = 0;

            //- VTK index
            virtual label mshIndex() const = 0;


        //- Set points based on vertices
        virtual tmp<pointField> getNodes(const List<vector>& verts) const = 0;

        //- Calculate shape functions
        virtual scalarList calcShape(const vector& pt) const = 0;

        //- Calculate the derivative of shape functions
        virtual scalarRectangularMatrix calcDShape(const vector& pt) const = 0;


    // Output functions

        //- Write to vtk
        virtual void vtkData
        (
            labelList& data,
            label& start,
            const labelList& labels
        ) const = 0;

        virtual const labelList& vtkConnectivity() const
        {
            if (!vtkConnectivity_.size())
            {
                vtkConnectivity_ =  vtkElementConnectivity(type_, order_);
            }
            return vtkConnectivity_;
        }

         //- Uniform points
        const List<point>& P() const
        {
            if (!P_.size())
            {
                createUniformGeometry();
            }
            return P_;
        }

        //- Uniform edges
        const List<labelPair>& E() const
        {
            if (!E_.size())
            {
                createUniformGeometry();
            }
            return E_;
        }

        //- Uniform geometry
        const List<labelList>& G() const
        {
            if (!G_.size())
            {
                createUniformGeometry();
            }
            return G_;
        }


    // Static member functions

        //- Return finite element pointer using ElementType
        static const finiteElement* getRefFiniteElement
        (
            const ElementType::Type et,
            const label order
        );

        //- Return finite element pointer using cell model name
        static const finiteElement* getRefFiniteElement
        (
            const word& cmt,
            const label order
        );

        static labelList matchNodes
        (
            const UList<vector>& localNodes,
            const List<vector>& globalNodes,
            const labelList& indices,
            const scalar tol = 1e-6
        );

        static labelList matchNodes
        (
            const UIndirectList<vector>& localNodes,
            const List<vector>& globalNodes,
            const labelList& indices,
            const scalar tol = 1e-6
        );
};

#define declareFE(Base, Order, Name)\
class Name : public Base            \
{                                   \
public:                             \
    TypeName(#Name)                 \
    Name() : Base(Order) {}         \
};

#define addFE(Name)                 \
defineTypeNameAndDebug(Name, 0);    \
addToRunTimeSelectionTable(finiteElement, Name, typeOrder);

#define addFEMap(Name, Table, ID)    \
addToRunTimeSelectionMap(finiteElement, Name, Table, ID);

class shellFiniteElement
{
public:

    typedef finiteElement::finiteElementTable finiteElementTable;

    //- Runtime type information
    TypeName("shellFiniteElement")

    // Declare runtime constructor selection tables

        //- Look up by type, given an order
        declareRunTimeSelectionTable
        (
            autoPtr,
            finiteElement,
            type,
            (
                const label order
            ),
            (order)
        );

        //- Look up by type, given an order
        declareRunTimeSelectionTable
        (
            autoPtr,
            finiteElement,
            typeOrder,
            (),
            ()
        );

    // Constructors

        //- Construct from order, size and element type
        shellFiniteElement();

    //- Destructor
    virtual ~shellFiniteElement();


    //- Stored shell elements
    static finiteElementTable finiteElements;

    // Static member functions

        //- Return finite element pointer using ElementType
        static const finiteElement* getRefFiniteElement
        (
            const ElementType::Type et,
            const label order
        );

        //- Return finite element pointer using cell model name
        static const finiteElement* getRefFiniteElement
        (
            const word& cmt,
            const label order
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "finiteElementI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
