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
    Foam::finiteElements::tetrahedron

Description
    Triangle finite element

SourceFiles
    tetrahedronFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef tetrahedronFiniteElement_H
#define tetrahedronFiniteElement_H

#include "segmentFiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class tetrahedron Declaration
\*---------------------------------------------------------------------------*/

class tetrahedron
:
    public FiniteElementBase<ElementType::TET>
{
    //- Reference to segment element
    const segment& seg_;

    //- Number of face nodes
    const label nFaceNodes_;

    //- Inverse transformation matrix
    scalarSquareMatrix invT_;
    labelList pivotIndices_;

    //- Create uniform geometry
    virtual void createUniformGeometry() const;

public:

    TypeName("tetrahedron")

    //- Constructor
    tetrahedron(const label order);

    //- Destructor
    virtual ~tetrahedron();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return (order_ + 1)*(order_ + 2)*(order_ + 3)/6;
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return order_ - 1;
            }

            //- Return the number of nodes for a face, neglecting edges
            inline label nFaceNodes(const label) const
            {
                label n = 0;
                for (label i = 2; i < order_; i++) n += order_ - i;
                return n;
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
        virtual tmp<pointField> getNodes(const List<vector>& verts) const;

        //- Calculate shape functions
        virtual scalarList calcShape(const vector& pt) const;

        //- Calculate the derivative of shape functions
        virtual scalarRectangularMatrix calcDShape(const vector& pt) const;


    // Output functions

        //- Write to vtk
        virtual void vtkData
        (
            labelList& data,
            label& start,
            const labelList& labels
        ) const;
};

declareFE(tetrahedron, 1, tet4);
declareFE(tetrahedron, 2, tet10);
declareFE(tetrahedron, 3, tet20);
declareFE(tetrahedron, 4, tet35);
declareFE(tetrahedron, 5, tet56);
declareFE(tetrahedron, 6, tet84);
declareFE(tetrahedron, 7, tet120);
declareFE(tetrahedron, 8, tet165);
declareFE(tetrahedron, 9, tet220);
declareFE(tetrahedron, 10, tet286);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
