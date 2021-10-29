/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technology
     \\/     M anipulation  |
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
    Foam::immersedTriangle

Description


SourceFiles
    immersedTriangle.C

\*---------------------------------------------------------------------------*/

#ifndef immersedTriangle_H
#define immersedTriangle_H

#include "immersedShape.H"
#include "FixedList.H"
#include "barycentric2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class immersedTriangle Declaration
\*---------------------------------------------------------------------------*/

class immersedTriangle
:
    public immersedShape
{
    // Private data

        //- List of points
        pointField vertices_;

        //- Return 2D points
        virtual void get2DPoints();


public:

    //- Runtime type information
    TypeName("triangle");


    // Constructors

        //- Construct from patch
        immersedTriangle
        (
            const polyMesh&,
            const immersedBoundaryObject&,
            const dictionary& dict
        );



    //- Destructor
    virtual ~immersedTriangle();


    // Member Functions

        //- Return bool list if list of points is inside the object
        virtual labelList calcInside(const pointField& points) const;

        //- Return if point is inside the object
        virtual bool inside(const point& pt) const;

        //- Return moment of inertia
        virtual diagTensor momentOfInertia() const
        {
            // Assumes that x1->x2 is along the x axis
            // This should have been handles in the constructor
            const vector& p1(vertices_[0]);
            const vector& p2(vertices_[1]);
            const vector& p3(vertices_[2]);

            vector p21(0.5*(p2 + p1));
            scalar b = mag(p2[xi_] - p1[xi_]);
            scalar h = mag(p3[yi_] - p21[yi_]);
            scalar a = mag(p3[xi_] - p1[xi_]);
            return diagTensor
            (
                b*pow3(h)/36.0,
                (pow3(b)*h - sqr(b)*h*a + b*h*sqr(a))/36.0,
                (b*pow3(h) + pow3(b)*h - sqr(b)*h*a + b*h*sqr(a))/36.0
            );
        }

        //- Does the centre of mass need to be read
        virtual bool needCentreOfMass() const
        {
            return false;
        }

        //- Does the orientation need to be read
        virtual bool needOrientation() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
