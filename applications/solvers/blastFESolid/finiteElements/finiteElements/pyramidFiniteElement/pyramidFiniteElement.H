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
    Foam::finiteElements::pyramid

Description
    Prism finite element

SourceFiles
    pyramidFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef pyramidFiniteElement_H
#define pyramidFiniteElement_H

#include "segmentFiniteElement.H"
#include "triangleFiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class pyramid Declaration
\*---------------------------------------------------------------------------*/

class pyramid
:
    public FiniteElement<ElementType::PYR>
{
    //- Create uniform geometry
    virtual void createUniformGeometry() const;

public:

    TypeName("pyramid")

    //- Constructor
    pyramid(const label order);

    //- Destructor
    virtual ~pyramid();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return 5;
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return 0;
            }

            //- Return the number of nodes for a face, neglecting edges
            inline label nFaceNodes(const label fi) const
            {
                return 0;
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

declareFE(pyramid, 1, pyramid5);
// declareFE(pyramid, 2, pyramid14);
// declareFE(pyramid, 3, pyramid30);
// declareFE(pyramid, 4, pyramid55);
// declareFE(pyramid, 5, pyramid91);
// declareFE(pyramid, 6, pyramid140);
// declareFE(pyramid, 7, pyramid204);
// declareFE(pyramid, 8, pyramid285);
// declareFE(pyramid, 9, pyramid385);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
