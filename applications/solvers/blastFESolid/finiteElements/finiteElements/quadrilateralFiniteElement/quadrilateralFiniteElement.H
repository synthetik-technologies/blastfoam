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
    Foam::finiteElements::quadrilateral

Description
    Quadrilateral finite element

SourceFiles
    quadrilateralFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef quadrilateralFiniteElement_H
#define quadrilateralFiniteElement_H

#include "segmentFiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class quadrilateral Declaration
\*---------------------------------------------------------------------------*/

class quadrilateral
:
    public FiniteElement<ElementType::QUAD>
{
    //- Reference to segment element
    const segment& seg_;

    //- Map from grid index to actual node index
    labelList map_;

    //- Create uniform geometry
    virtual void createUniformGeometry() const;


public:

    TypeName("quadrilateral")

    //- Constructor
    quadrilateral(const label order);

    //- Destructor
    virtual ~quadrilateral();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return (order_ + 1)*(order_ + 1);
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return order_ - 1;
            }

            //- Return the number of nodes for a face, neglecting edges
            inline label nFaceNodes(const label) const
            {
                return (order_ - 1)*(order_ - 1);
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

declareFE(quadrilateral, 1, quad4);
declareFE(quadrilateral, 2, quad9);
declareFE(quadrilateral, 3, quad16);
declareFE(quadrilateral, 4, quad25);
declareFE(quadrilateral, 5, quad36);
declareFE(quadrilateral, 6, quad49);
declareFE(quadrilateral, 7, quad64);
declareFE(quadrilateral, 8, quad81);
declareFE(quadrilateral, 9, quad100);
declareFE(quadrilateral, 10, quad121);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
