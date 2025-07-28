/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    functionObjects::centrifugalBodyForce

Description
    FunctionObject applys centrifugal body force

Author
    Zeljko Tukovic, FSB.  All rights reserved

SourceFiles
    centrifugalBodyForce.C

\*---------------------------------------------------------------------------*/

#ifndef centrifugalBodyForce_H
#define centrifugalBodyForce_H

#include "fvMeshFunctionObject.H"
#include "dictionary.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class centrifugalBodyForce Declaration
\*---------------------------------------------------------------------------*/

class centrifugalBodyForce
:
    public fvMeshFunctionObject
{
    // Private data

        //- Angular velocity
        dimensionedVector angularVelocity_;


public:

    //- Runtime type information
    TypeName("centrifugalBodyForce");

    // Constructors

        //- Construct from components
        centrifugalBodyForce
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
