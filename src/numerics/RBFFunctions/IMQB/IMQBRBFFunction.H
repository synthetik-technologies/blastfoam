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
    IMQB

Description
    IMQB radial basis function

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    IMQB.C

\*---------------------------------------------------------------------------*/

#ifndef IMQB_H
#define IMQB_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFFunctions
{

/*---------------------------------------------------------------------------*\
                              Class IMQB Declaration
\*---------------------------------------------------------------------------*/

class IMQB
:
    public RBFFunction
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IMQB(const IMQB&);

        //- Disallow default bitwise assignment
        void operator=(const IMQB&);


public:

    //- Runtime type information
    TypeName("IMQB");

    // Constructors

        //- Construct given radius
        IMQB(const scalar radius);

        //- Construct from dictionary
        IMQB(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new IMQB(this->radius_));
        }


    // Destructor

        virtual ~IMQB();


    // Member Functions

        //- Return weights given points
        virtual scalar evaluate(const scalar dist) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
