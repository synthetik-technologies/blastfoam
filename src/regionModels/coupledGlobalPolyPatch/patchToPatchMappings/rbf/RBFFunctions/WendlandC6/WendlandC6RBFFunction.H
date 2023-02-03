/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | Copyright David Blom
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derived work of foam-extend.

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
    Foam::RBFFunctions::WendlandC6

Description

Author
    David Blom, TU Delft. All rights reserved.

SourceFiles
    WendlandC6RBFFunction.C

\*---------------------------------------------------------------------------*/

#ifndef WendlandC6RBFFunction_H
#define WendlandC6RBFFunction_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFFunctions
{

/*---------------------------------------------------------------------------*\
                              Class WendlandC24 Declaration
\*---------------------------------------------------------------------------*/

class WendlandC6
:
    public RBFFunction
{
    // Private data

        //- Radius
        scalar radius_;


public:

    //- Runtime type information
    TypeName("WendlandC6");

    // Constructors

        //- Construct from dictionary
        WendlandC6(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new WendlandC6(*this));
        }


    // Destructor

        virtual ~WendlandC6();


    // Member Functions

        inline scalar evaluate(const scalar dist) const
        {
            const scalar p(dist/radius_);
            return
                p > 1
              ? 0.0
              : sqr(pow4(1.0 - p))
               *(32.0*pow3(p) + 25.0*sqr(p) + 8.0*p + 1.0);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
