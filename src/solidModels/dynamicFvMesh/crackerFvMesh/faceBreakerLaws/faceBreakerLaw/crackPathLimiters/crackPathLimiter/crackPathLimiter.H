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
    crackPathLimiter

Description
    Law that limits where cracks are allowed to propagate.

    Returns surfaceScalarField with 1's and 0's;
        1 means face is allowed to break
        0 means face is not allowed to break

Authors
    Philip Cardif, UCD. All rights reserved.

SourceFiles
    crackPathLimiter.C
    newCrackPathLimiter.C

\*---------------------------------------------------------------------------*/

#ifndef crackPathLimiter_H
#define crackPathLimiter_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class crackPathLimiter Declaration
\*---------------------------------------------------------------------------*/

class crackPathLimiter
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        crackPathLimiter(const crackPathLimiter&);

        //- Disallow default bitwise assignment
        void operator=(const crackPathLimiter&);


protected:

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

public:

    //- Runtime type information
    TypeName("crackPathLimiter")


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            crackPathLimiter,
            dictionary,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected meshFailure model
        static autoPtr<crackPathLimiter> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        crackPathLimiter
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~crackPathLimiter()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        // Return list of faces allowed to break
        virtual const surfaceScalarField& facesAllowedToBreak() const = 0;

        //- Delete demand driven data
        virtual void clearOut() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
