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
    Foam::finiteElements::hexahedron

Description
    Hexahedron finite element

SourceFiles
    hexahedronFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef hexahedronFiniteElement_H
#define hexahedronFiniteElement_H

#include "segmentFiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class hexahedron Declaration
\*---------------------------------------------------------------------------*/

class hexahedron
:
    public FiniteElement<ElementType::HEX>
{
    //- Reference to segment element
    const segment& seg_;

    //- Map from grid index to actual node index
    labelList map_;

    //- Create uniform geometry
    virtual void createUniformGeometry() const;

public:

    TypeName("hexahedron")

    //- Constructor
    hexahedron(const label order);

    //- Destructor
    virtual ~hexahedron();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return (order_ + 1)*(order_ + 1)*(order_ + 1);
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
                return (order_ - 1)*(order_ - 1)*(order_ - 1);
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

declareFE(hexahedron, 1, hex8);
declareFE(hexahedron, 2, hex27);
declareFE(hexahedron, 3, hex64);
declareFE(hexahedron, 4, hex125);
declareFE(hexahedron, 5, hex216);
declareFE(hexahedron, 6, hex343);
declareFE(hexahedron, 7, hex512);
declareFE(hexahedron, 8, hex729);
declareFE(hexahedron, 9, hex1000);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
