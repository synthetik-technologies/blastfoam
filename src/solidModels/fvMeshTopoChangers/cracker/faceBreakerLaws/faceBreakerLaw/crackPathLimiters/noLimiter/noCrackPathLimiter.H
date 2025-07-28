/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::crackPathLimiters::none

Description
    Limit cracks to propagate inside specified bounding boxes.

    Usage:
    verbatim
        crackPathLimiterCoeffs
        {
            // Cracks may only propagate along faces within these boxes
            // This allows the crack path to be limited.
            crackPathLimiter    nonees;
            nonees
            (
                (0 -1e-5 0) (0.05 1e-5 1)
                (0.05 1e-5 0) (0.05 2e-5 1)
            );
        }
    endverbatim

Author
    Philip Cardif, UCD.  All rights reserved.

SourceFiles
    noneCrackPathLimiter.C

\*---------------------------------------------------------------------------*/

#ifndef noneCrackPatchLimiter_H
#define noneCrackPatchLimiter_H

#include "crackPathLimiter.H"
#include "boundBox.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace crackPathLimiters
{

/*---------------------------------------------------------------------------*\
                         Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public crackPathLimiter
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        none(const none&);

        //- Disallow default bitwise assignment
        void operator=(const none&);


public:

    //- Runtime type information
    TypeName("none");

    // Static data members


    // Constructors

        //- Construct from dictionary
        none
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~none();


    // Member Functions

        // Return list of faces allowed to break
        virtual tmp<surfaceScalarField> facesAllowedToBreak() const;

        //- Delete demand driven data
        virtual void clearOut();

        //- Write
        virtual bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace crackPathLimiters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
