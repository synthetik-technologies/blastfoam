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
    Foam::finiteElements::prism

Description
    Prism finite element

SourceFiles
    prismFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef prismFiniteElement_H
#define prismFiniteElement_H

#include "segmentFiniteElement.H"
#include "triangleFiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class prism Declaration
\*---------------------------------------------------------------------------*/

class prism
:
    public FiniteElementBase<ElementType::PRISM>
{
    //- Reference to segment element
    const segment& seg_;

    //- Reference to triangle element
    const triangle& tri_;

    //- Create uniform geometry
    virtual void createUniformGeometry() const;

public:

    TypeName("prism")

    //- Constructor
    prism(const label order);

    //- Destructor
    virtual ~prism();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return (order_ + 1)*(order_ + 1)*(order_ + 2)/2.0;
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return order_ - 1;
            }

            //- Return the number of nodes for a face, neglecting edges
            inline label nFaceNodes(const label fi) const
            {
                return
                    fi < 2
                  ? tri_.nFaceNodes(fi)
                  : (order_ - 1)*(order_ - 1);
            }

            //- Return the number of internal nodes
            inline label nInternalNodes() const
            {
                return tri_.nFaceNodes(0)*(order_ - 1);
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

declareFE(prism, 1, wedge6);
declareFE(prism, 2, wedge18);
declareFE(prism, 3, wedge40);
declareFE(prism, 4, wedge75);
declareFE(prism, 5, wedge126);
declareFE(prism, 6, wedge196);
declareFE(prism, 7, wedge288);
declareFE(prism, 8, wedge405);
declareFE(prism, 9, wedge550);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
