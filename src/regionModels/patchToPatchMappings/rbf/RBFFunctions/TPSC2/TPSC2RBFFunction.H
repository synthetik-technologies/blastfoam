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
    TPSC2

Description
    TPSC2 radial basis function

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    TPSC2RBFFunction.C

\*---------------------------------------------------------------------------*/

#ifndef TPSC2RBFFunction_H
#define TPSC2RBFFunction_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFFunctions
{

/*---------------------------------------------------------------------------*\
                              Class TPSC2 Declaration
\*---------------------------------------------------------------------------*/

class TPSC2
:
    public RBFFunction
{
    // Private data

        //- Inverse support radius
        scalar invRadius_;


public:

    //- Runtime type information
    TypeName("TPSC2");

    // Constructors

        //- Construct from dictionary
        TPSC2(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new TPSC2(*this));
        }


    // Destructor

        virtual ~TPSC2();


    // Member Functions

        //- Return weights given points
        inline scalar evaluate(const scalar dist) const
        {
            const scalar p = dist*invRadius_;
            return
                1.0
              - 30.0*sqr(p)
              - 10.0*pow3(p)
              + 45.0*pow4(p)
              - 6.0*pow5(p)
              - 60.0*pow3(p)*log(max(p, small));
        }

        //- Return the cutoff radius
        virtual scalar cutOffRadius() const
        {
            return 1.0/invRadius_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
