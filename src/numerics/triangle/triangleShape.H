/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::triangle

Description
    A triangle primitive used to calculate face normals and swept volumes.

SourceFiles
    triangleShapeI.H

\*---------------------------------------------------------------------------*/

#ifndef triangleShape_H
#define triangleShape_H

#include "triangle.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
class Istream;

// Forward declaration of friend functions and operators

template<class Point, class PointRef> class triangleShape;

template<class Point, class PointRef>
inline Istream& operator>>
(
    Istream&,
    triangleShape<Point, PointRef>&
);


/*---------------------------------------------------------------------------*\
                           class triangle Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class triangleShape
:
    public triangle<Point, PointRef>
{

public:

    // Constructors

        //- Construct from three points
        inline triangleShape(const Point& a, const Point& b, const Point& c);

        //- Construct from Istream
        inline triangleShape(Istream&);


    // Member Functions

        // Properties

            //- Return shape function gradient
            inline vector gradNi(const label i) const;

            //- Evaluate shape function at a given point
            inline scalar Ni(const label i, const Point& p) const;


    // IOstream operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            triangleShape&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triangleShapeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
