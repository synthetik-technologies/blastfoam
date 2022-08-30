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
    epsilonPEqCellRemovalLaw

Description
    Cell removal law which selects cells with an equivalent plastic strain
    greater than the specified critical equivalent value.

Author
    Philip Cardif, UCD/UT. All rights reserved.

SourceFiles
    epsilonPEqCellRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonPEqCellRemovalLaw_H
#define epsilonPEqCellRemovalLaw_H

#include "cellRemovalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class epsilonPEqCellRemovalLaw Declaration
\*---------------------------------------------------------------------------*/

class epsilonPEqCellRemovalLaw
:
    public cellRemovalLaw
{
    // Private Data Member

        // Critical value of epsilonPEq
        const scalar epsilonPEqCrit_;

        // Name of the equivalent plastic strain field
        const word epsilonPEqName_;

        // Index of patch to put exposed faces in
        const scalar patchID_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        epsilonPEqCellRemovalLaw(const epsilonPEqCellRemovalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const epsilonPEqCellRemovalLaw&);


public:

    //- Runtime type information
    TypeName("epsilonPEq");

    // Static data members


    // Constructors

        //- Construct from dictionary
        epsilonPEqCellRemovalLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~epsilonPEqCellRemovalLaw();


    // Member Functions

        // Return list of cell IDs to be removed
        virtual tmp<labelField> cellsToRemove();

        // Return patch index for the faces to be exposed
        virtual label exposedFacesPatchID();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
