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
    Foam::finiteElements::triangle

Description
    Triangle finite element

SourceFiles
    triangleFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef triangleFiniteElement_H
#define triangleFiniteElement_H

#include "segmentFiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class triangle Declaration
\*---------------------------------------------------------------------------*/

class triangle
:
    public FiniteElement<ElementType::TRI>
{
    //- Reference to segment element
    const segment& seg_;

    const label nFaceNodes_;

    //- Inverse transformation matrix
    scalarSquareMatrix invT_;

    //- Create uniform geometry
    virtual void createUniformGeometry() const;


public:

    TypeName("triangle")

    //- Constructor
    triangle(const label order);

    //- Destructor
    virtual ~triangle();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return (order_ + 1)*(order_ + 2)/2;
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return order_ - 1;
            }

            //- Return the number of nodes for a face, neglecting edges
            inline label nFaceNodes(const label) const
            {
                return nFaceNodes_;
            }

            //- Return the number of internal nodes
            inline label nInternalNodes() const
            {
                return 0;
            }

            //- VTK index
            virtual label vtkIndex() const;

            //- VTK index
            virtual label mshIndex() const;


        //- Set points based on vertices
        virtual tmp<pointField> getNodes
        (
            const vector& p0,
            const vector& p1,
            const vector& p2
        ) const;

        //- Set points based on vertices
        virtual tmp<pointField> getNodes(const List<vector>& verts) const;

        //- Calculate shape functions
        virtual scalarList calcShape(const vector& pt) const;

        //- Calculate the derivative of shape functions
        virtual scalarRectangularMatrix calcDShape(const vector& pt) const;

        //- Calculate the shape function and its derivatives
        virtual void calcDShape
        (
            const vector& pt,
            List<scalar>& shape,
            scalarRectangularMatrix& dshape
        ) const;


    // Output functions

        //- Write to vtk
        virtual void vtkData
        (
            labelList& data,
            label& start,
            const labelList& labels
        ) const;


    // Static member functions

        //- Calculate the number of nodes for a face, neglecting edges
        static label calcNFaceNodes(const label o);
};

declareFE(triangle, 1, tri3);
declareFE(triangle, 2, tri6);
declareFE(triangle, 3, tri10);
declareFE(triangle, 4, tri15);
declareFE(triangle, 5, tri21);
declareFE(triangle, 6, tri28);
declareFE(triangle, 7, tri36);
declareFE(triangle, 8, tri45);
declareFE(triangle, 9, tri55);
declareFE(triangle, 10, tri66);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
