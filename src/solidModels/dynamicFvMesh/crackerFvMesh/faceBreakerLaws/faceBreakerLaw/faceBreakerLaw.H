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
    faceBreakerLaw

Description
    Law that selects which internal faces should be cracked/broken.

Authors
    Philip Cardif, UCD. All rights reserved.

SourceFiles
    faceBreakerLaw.C
    newCellRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef faceBreakerLaw_H
#define faceBreakerLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faceBreakerLaw Declaration
\*---------------------------------------------------------------------------*/

class faceBreakerLaw
{
    // Private data

        //- Name
        const word name_;

        //- Reference to the field
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        faceBreakerLaw(const faceBreakerLaw&);

        //- Disallow default bitwise assignment
        void operator=(const faceBreakerLaw&);


protected:

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

public:

    //- Runtime type information
    TypeName("faceBreakerLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceBreakerLaw,
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
        static autoPtr<faceBreakerLaw> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        faceBreakerLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~faceBreakerLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        // Return list of internal faces to break
        virtual const labelList& facesToBreak() const = 0;

        // Return list of coupled faces to break
        virtual const labelList& coupledFacesToBreak() const = 0;

        //- Delete demand driven data
        virtual void clearOut() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
