/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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
    Foam::element

Description
    Element for use with feMesh

SourceFiles
    element.C

\*---------------------------------------------------------------------------*/
#ifndef element_H
#define element_H

#include "finiteElement.H"
#include "GeoType.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class primitiveMesh;
class element;

Istream& operator>>(Istream& is, element& elem);
Ostream& operator>>(Ostream& is, const element& elem);

/*---------------------------------------------------------------------------*\
                        Class element Declaration
\*---------------------------------------------------------------------------*/

class element
:
    public List<label>
{
public:

    //- Enumeration of read type
    enum ReadType
    {
        VTK,
        MSH
    };

protected:

    //- Pointer to the reference finite element
    const finiteElement* elem_;

    //- Integration rule for this element
    const integrationRule* ir_;

    //- Geometric type
    GeoType::Type geoType_;


public:

    // Constructors

        //- Null constructor
        element();

        //- Construct from mesh, index, order, and geometric type
        element
        (
            const primitiveMesh& mesh,
            const label index,
            const label order,
            const GeoType::Type type
        );

        //- Construct from a read element
//         element(const labelList& elem, const ReadType rt);


    //- Destructor
    ~element();


    // Member functions

        // Access

            //- Return geometric type
            inline GeoType::Type geoType() const;

            //- Return the element type
            inline ElementType::Type type() const;

            //- Return the number of dimensions
            inline label nDims() const;

            //- Return element order
            inline label order() const;

            //- Return the integration rule
            inline const integrationRule& ir() const;

            //- Return number of vertices
            inline label nVertices() const;

            //- Return number of edges
            inline label nEdges() const;

            //- Return number of faces
            inline label nFaces() const;

            //- Return number of nodes
            inline label nNodes() const;

            //- Return number of edge nodes (neglecting ends)
            inline label nEdgeNodes() const;

            //- Return number of face nodes (neglecting edges)
            inline label nFaceNodes(const label fi) const;

            //- Return number of internal nodes
            inline label nInternalNodes() const;

            //- Return reference nodes
            inline const pointField& nodes() const;

            //- Return the reference finite element
            inline const finiteElement& elem() const;

            // Return edges mapped from this element
            inline edgeList edges() const;

            // Return faces mapped from this element
            inline faceList faces() const;

            //- Return all nodes in physical space
            inline tmp<pointField> getNodes(const List<vector>& pts) const;


            //- Set the integtaion rule and update save data
            inline void setIntegrationRule(const integrationRule& ir);

            //- Set the reference finiteElement
            inline void setFiniteElement(const finiteElement& fe);

            //- Set the reference finiteElement
            inline void setGeoType(const GeoType::Type& gt);

        //- Construct from mesh, index, order, and geometric type
        void initialize
        (
            const primitiveMesh& mesh,
            const label index,
            const label order,
            const GeoType::Type type
        );

        //- Construct from a read element
//         void initialize(const labelList& elem, const ReadType rt);

        //- Calculate the shape functions at an integration point
        scalarList calcShape(const vector& pt) const;

        //- Calculate the derivates of shape functions at an
        //  integration point
        scalarRectangularMatrix calcDShape(const vector& pt) const;

        //- Calculate the Jacobian at an integration point
        tensor calcJ(const vector& ip, const pointField& meshNodes) const;


        //- Calculate the inverse Jacobian at an integration point
        tensor calcInvJ(const vector& ip, const pointField& meshNodes) const;

        //- Calculate the inverse Jacobian given the current Jacobian
        tensor calcInvJ(const tensor& J) const;

        //- Calculate the weight at an integration point
        scalar calcW(const vector& ip, const pointField& meshNodes) const;

        //- Calculate the weight given the current Jacobian
        scalar calcW(const tensor& J) const;


    // IO operators

        friend Istream& operator>>(Istream& is, element& elem);
        friend Ostream& operator>>(Ostream& is, const element& elem);
};

typedef List<element> elementList;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "elementI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
