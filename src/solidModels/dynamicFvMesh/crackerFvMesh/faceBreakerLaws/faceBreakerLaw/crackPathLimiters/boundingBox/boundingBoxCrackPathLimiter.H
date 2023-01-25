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
    Foam::crackPathLimiters::boundingBox

Description
    Limit cracks to propagate inside specified bounding boxes.

    Usage:
    verbatim
        crackPathLimiterCoeffs
        {
            // Cracks may only propagate along faces within these boxes
            // This allows the crack path to be limited.
            crackPathLimiter    boundingBoxes;
            boundingBoxes
            (
                (0 -1e-5 0) (0.05 1e-5 1)
                (0.05 1e-5 0) (0.05 2e-5 1)
            );
        }
    endverbatim

Author
    Philip Cardif, UCD.  All rights reserved.

SourceFiles
    boundingBoxCrackPathLimiter.C

\*---------------------------------------------------------------------------*/

#ifndef boundingBoxCrackPatchLimiter_H
#define boundingBoxCrackPatchLimiter_H

#include "crackPathLimiter.H"
#include "boundBox.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace crackPathLimiters
{
/*---------------------------------------------------------------------------*\
                         Class boundingBox Declaration
\*---------------------------------------------------------------------------*/

class boundingBox
:
    public crackPathLimiter
{
    // Private data

        // faces allowed to break: 1 is allowed, 0 is not allowed
        mutable surfaceScalarField* facesAllowedToBreakPtr_;

        // List of bounding boxes
        List<boundBox> boundingBoxes_;


    // Private Member Functions

        // Calculate faces allowed to break
        void calcFacesAllowedToBreak() const;

        //- Disallow default bitwise copy construct
        boundingBox(const boundingBox&);

        //- Disallow default bitwise assignment
        void operator=(const boundingBox&);


public:

    //- Runtime type information
    TypeName("boundingBoxes");

    // Static data members


    // Constructors

        //- Construct from dictionary
        boundingBox
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~boundingBox();


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
