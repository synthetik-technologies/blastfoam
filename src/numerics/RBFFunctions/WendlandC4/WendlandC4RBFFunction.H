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
    Foam::RBFFunctions::WendlandC4

Description

Author
    David Blom, TU Delft. All rights reserved.

SourceFiles
    WendlandC4RBFFunction.C

\*---------------------------------------------------------------------------*/

#ifndef WendlandC4RBFFunction_H
#define WendlandC4RBFFunction_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFFunctions
{

/*---------------------------------------------------------------------------*\
                              Class WendlandC24 Declaration
\*---------------------------------------------------------------------------*/

class WendlandC4
:
    public RBFFunction
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WendlandC4(const WendlandC4&);

        //- Disallow default bitwise assignment
        void operator=(const WendlandC4&);


public:

    //- Runtime type information
    TypeName("WendlandC4");

    // Constructors

        //- Construct given radius
        WendlandC4(const scalar radius);

        //- Construct from dictionary
        WendlandC4(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new WendlandC4(this->radius_));
        }


    // Destructor

        virtual ~WendlandC4();


    // Member Functions

        virtual scalar evaluate(const scalar dist) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
