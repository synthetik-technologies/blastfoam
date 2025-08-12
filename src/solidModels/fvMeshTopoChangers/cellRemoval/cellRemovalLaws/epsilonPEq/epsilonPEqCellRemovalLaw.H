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
    Foam::cellRemovalLaws::epsilonPEq

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
namespace cellRemovalLaws
{

/*---------------------------------------------------------------------------*\
                         Class epsilonPEqCell Declaration
\*---------------------------------------------------------------------------*/

class epsilonPEq
:
    public cellRemovalLaw
{
    // Private Data Member

        // Critical value of epsilonPEq
        const scalar epsilonPEqCrit_;

        // Name of the equivalent plastic strain field
        const word epsilonPEqName_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        epsilonPEq(const epsilonPEq&);

        //- Disallow default bitwise assignment
        void operator=(const epsilonPEq&);


public:

    //- Runtime type information
    TypeName("epsilonPEq");

    // Static data members


    // Constructors

        //- Construct from dictionary
        epsilonPEq
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~epsilonPEq();


    // Member Functions

        // Return list of cell IDs to be removed
        virtual labelHashSet cellsToRemove();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellRemovalLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
