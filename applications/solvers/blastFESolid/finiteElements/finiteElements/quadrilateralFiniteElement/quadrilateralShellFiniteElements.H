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
    Foam::finiteElements::quadrilateral4Shell

Description
    Quadrilateral finite element

SourceFiles
    quadrilateral4ShellFiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef quadrilateral4ShellFiniteElement_H
#define quadrilateral4ShellFiniteElement_H

#include "finiteElement.H"
#include "FiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class quadrilateral4Shell Declaration
\*---------------------------------------------------------------------------*/

class quadrilateral4Shell
:
    public FiniteElement<ElementType::QUAD>
{
    //- Create uniform geometry
    virtual void createUniformGeometry() const
    {
        NotImplemented;
    }


public:

    //- Constructor
    quadrilateral4Shell();

    //- Destructor
    virtual ~quadrilateral4Shell();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return 4;
            }

            //- Return the number of edge nodes (neglecting start/end)
            inline label nEdgeNodes() const
            {
                return 0;
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
            inline label vtkIndex() const
            {
                return 9;
            }

            //- VTK index
            inline label mshIndex() const
            {
                return 3;
            }


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


/*---------------------------------------------------------------------------*\
                        Class quadrilateral4Shell1 Declaration
\*---------------------------------------------------------------------------*/

class quadrilateral4Shell1
:
    public quadrilateral4Shell
{
public:

    TypeName("quad1Shell1")

    //- Constructor
    quadrilateral4Shell1();
};


/*---------------------------------------------------------------------------*\
                        Class quadrilateral4Shell2 Declaration
\*---------------------------------------------------------------------------*/

class quadrilateral4Shell2
:
    public quadrilateral4Shell
{
public:

    TypeName("quad1Shell2")

    //- Constructor
    quadrilateral4Shell2();
};


/*---------------------------------------------------------------------------*\
                        Class quadrilateral4Shell3 Declaration
\*---------------------------------------------------------------------------*/

class quadrilateral4Shell3
:
    public quadrilateral4Shell
{
public:

    TypeName("quad1Shell3")

    //- Constructor
    quadrilateral4Shell3();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
