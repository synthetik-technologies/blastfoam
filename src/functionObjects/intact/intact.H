/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::intact

Description
    Write out a field corresponding to intact patches, used with burst patches

Usage
    \table
        Property          | Description               | Required | Default
        type              | type name: intact         | yes      |


See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    intact.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_intact_H
#define functionObjects_intact_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class intact Declaration
\*---------------------------------------------------------------------------*/

class intact
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("intact");


    // Constructors

        //- Construct from Time and dictionary
        intact
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        intact(const intact&) = delete;


    //- Destructor
    virtual ~intact();


    // Member Functions

        //- Calculate the field averages
        virtual bool execute()
        {
            return true;
        }

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const intact&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
