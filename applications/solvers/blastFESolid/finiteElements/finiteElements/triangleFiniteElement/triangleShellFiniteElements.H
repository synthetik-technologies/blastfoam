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
    Foam::finiteElements::triangle3Shell
    Foam::finiteElements::triangle3Shell1
    Foam::finiteElements::triangle3Shell2
    Foam::finiteElements::triangle3Shell3

Description
    Three node triangle shell finite elements

SourceFiles
    triangle3ShellFiniteElements.C

\*---------------------------------------------------------------------------*/

#ifndef triangle3ShellFiniteElements_H
#define triangle3ShellFiniteElements_H

#include "finiteElement.H"
#include "FiniteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{

/*---------------------------------------------------------------------------*\
                        Class triangle3Shell Declaration
\*---------------------------------------------------------------------------*/

class triangle3Shell
:
    public FiniteElement<ElementType::TRI>,
    public shellFiniteElement
{
public:

    //- Stored integration rules
    static List<integrationRule> integrationRules;

protected:

    //- Create uniform geometry
    virtual void createUniformGeometry() const
    {
        NotImplemented;
    }

public:

    //- Constructor
    triangle3Shell();

    //- Destructor
    virtual ~triangle3Shell();


    // Member functions

        // Access

            //- Return the number of nodes
            inline label nNodes() const
            {
                return 3;
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
                return 5;
            }

            //- VTK index
            inline label mshIndex() const
            {
                return 1;
            }

        //- Return integration rule for a given order
        virtual const integrationRule& ir
        (
            const label orderi,
            const label orderj,
            const label orderk
        ) const;

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
};


/*---------------------------------------------------------------------------*\
                        Class triangle3Shell1 Declaration
\*---------------------------------------------------------------------------*/

class triangle3Shell1
:
    public triangle3Shell
{
public:

    TypeName("tri1Shell1")

    //- Constructor
    triangle3Shell1();
};


/*---------------------------------------------------------------------------*\
                        Class triangle3Shell2 Declaration
\*---------------------------------------------------------------------------*/

class triangle3Shell2
:
    public triangle3Shell
{
public:

    TypeName("tri1Shell2")

    //- Constructor
    triangle3Shell2();
};


/*---------------------------------------------------------------------------*\
                        Class triangle3Shell3 Declaration
\*---------------------------------------------------------------------------*/

class triangle3Shell3
:
    public triangle3Shell
{
public:

    TypeName("tri1Shell3")

    //- Constructor
    triangle3Shell3();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace finiteElements
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
