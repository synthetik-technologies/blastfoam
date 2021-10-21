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


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WendlandC6(const WendlandC6&);

        //- Disallow default bitwise assignment
        void operator=(const WendlandC6&);


public:

    //- Runtime type information
    TypeName("WendlandC6");

    // Constructors

        //- Construct given radius
        WendlandC6(const scalar radius);

        //- Construct from dictionary
        WendlandC6(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new WendlandC6(this->radius_));
        }


    // Destructor

        virtual ~WendlandC6();


    // Member Functions

        virtual scalar evaluate(const scalar dist) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
