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
    Foam::cellRemovalLaws::stress

Description
    Cell removal law which selects cells with an stress
    greater than the specified critical value.

SourceFiles
    stressCellRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef stressCellRemovalLaw_H
#define stressCellRemovalLaw_H

#include "cellRemovalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellRemovalLaws
{

/*---------------------------------------------------------------------------*\
                         Class stress Declaration
\*---------------------------------------------------------------------------*/

class stress
:
    public cellRemovalLaw
{
    // Private Data Member

        //- use von Mises stress measure
        bool vonMises_;

        // Critical value of sigma
        const scalar stressCrit_;

        // Name of the equivalent plastic strain field
        const word sigmaName_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stress(const stress&);

        //- Disallow default bitwise assignment
        void operator=(const stress&);


public:

    //- Runtime type information
    TypeName("stress");

    // Constructors

        //- Construct from dictionary
        stress
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~stress();


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
