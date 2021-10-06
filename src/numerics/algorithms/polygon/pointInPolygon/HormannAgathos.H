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
    HormannAgathos

Description
    Implements the point in polygon problem using the winding number technique
    presented in the paper:
        "The point in polygon problem for arbitrary polygons",
        Kai Hormann, Alexander Agathos, 2001

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

SourceFiles
    HormannAgathosI.H
    HormannAgathos.C
    HormannAgathosIO.C

\*---------------------------------------------------------------------------*/

#ifndef HormannAgathos_H
#define HormannAgathos_H

#include "List.H"
#include "point2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class HormannAgathos Declaration
\*---------------------------------------------------------------------------*/

class HormannAgathos
{
    // Private data

        //- 2D coordinate of polygon vertices
        //  We keep the same notation as in the paper
        List<point2D> P_;

        //- 2D distance tolerance factor for in/out tests
        scalar distTol_;

        //- 2D distance epsilon for in/out tests
        scalar epsilon_;


    // Private Member Functions

        //- Compute 2D distance epsilon based on a tolerance factor
        void evaluateEpsilon();

        //- Comparison of two scalar within a tolerance
        inline bool equalWithTol
        (
            const scalar& a,
            const scalar& b
        ) const;

        inline bool greaterWithTol
        (
            const scalar& a,
            const scalar& b
        ) const;

        inline bool smallerWithTol
        (
            const scalar& a,
            const scalar& b
        ) const;

        inline bool greaterOrEqualWithTol
        (
            const scalar& a,
            const scalar& b
        ) const;

        inline bool smallerOrEqualWithTol
        (
            const scalar& a,
            const scalar& b
        ) const;


public:

    // Public typedefs
    enum inOutClassification
    {
        POINT_OUTSIDE,
        POINT_INSIDE,
        POINT_ON_VERTEX,
        POINT_ON_EDGE
    };

    // Constructors

        //- Construct from components
        HormannAgathos
        (
            const List<point2D>& P,
            const scalar& distTol
        );


    // Destructor - default


    // Member Functions

        //- Executa classification of points
        inOutClassification evaluate(const point2D& R) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HormannAgathosI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
