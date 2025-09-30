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
    Foam::finiteElements::segment

Description
    Segment finite element

SourceFiles
    pointFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef segmentFiniteElement_H
#define segmentFiniteElement_H

#include "FiniteElementBase.H"
#include "CmptList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class segment Declaration
\*---------------------------------------------------------------------------*/

class segment
:
    public FiniteElementBase<ElementType::SEG>
{
    //- Save weights
    List<scalar> w_;

    //- Create uniform geometry
    virtual void createUniformGeometry() const;


public:

    TypeName("segment")

    //- Constructor
    segment(const label order);

    //- Destructor
    virtual ~segment();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return order_ + 1;
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return order_ - 1;
            }

            //- Return the number of nodes for a face, neglecting edges
            inline label nFaceNodes(const label) const
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

        //- Calculate shape functions, no allocation
        void calcShape(const scalar x, scalar* shape) const;

        //- Calculate shape functions
        scalarList calcShape(const scalar x) const;

        //- Calculate shape functions
        virtual scalarList calcShape(const vector& pt) const;

        //- Calculate the derivative of shape functions (no allocation)
        void calcDShape(const scalar x, scalar*) const;

        //- Calculate the derivative of shape functions
        scalarList calcDShape(const scalar x) const;

        //- Calculate the shape function and derivatives
        void calcDShape
        (
            const scalar x,
            List<scalar>& shape,
            List<scalar>& dshape
        ) const;

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

declareFE(segment, 1, seg2);
declareFE(segment, 2, seg3);
declareFE(segment, 3, seg4);
declareFE(segment, 4, seg5);
declareFE(segment, 5, seg6);
declareFE(segment, 6, seg7);
declareFE(segment, 7, seg8);
declareFE(segment, 8, seg9);
declareFE(segment, 9, seg10);
declareFE(segment, 10, seg11);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
