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
    Foam::RBFFunctions::WendlandC0

Description

Author
    David Blom, TU Delft. All rights reserved.

SourceFiles
    WendlandC0RBFFunction.C

\*---------------------------------------------------------------------------*/

#ifndef WendlandC0RBFFunction_H
#define WendlandC0RBFFunction_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFFunctions
{

/*---------------------------------------------------------------------------*\
                              Class WendlandC0 Declaration
\*---------------------------------------------------------------------------*/

class WendlandC0
:
    public RBFFunction
{
    // Private data

        //- Radius
        scalar radius_;


public:

    //- Runtime type information
    TypeName("WendlandC0");

    // Constructors

        //- Construct from dictionary
        WendlandC0(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new WendlandC0(*this));
        }


    // Destructor

        virtual ~WendlandC0();


    // Member Functions

        //- Return weights given points
        inline scalar evaluate(const scalar dist) const
        {
            return sqr(1.0 - dist/radius_);
        }

        //- Return the cutoff radius
        virtual scalar cutOffRadius() const
        {
            return radius_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
