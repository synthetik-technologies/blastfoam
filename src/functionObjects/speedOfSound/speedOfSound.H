/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate speed of sound with blastFoam thermo
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
    Foam::functionObjects::speedOfSound

Description
    Calculates and writes the dynamic pressure. Includes option to save the field
    to the database.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    speedOfSound.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_speedOfSound_H
#define functionObjects_speedOfSound_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class speedOfSound Declaration
\*---------------------------------------------------------------------------*/

class speedOfSound
:
    public fvMeshFunctionObject
{
    // Private Member data

        //- Name of phase
        word phaseName_;

        //- Name of Result field
        word resultName_;

        //- Is the field stored in the database
        Switch store_;


public:

    //- Runtime type information
    TypeName("speedOfSound");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        speedOfSound
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        speedOfSound(const speedOfSound&) = delete;


    //- Destructor
    virtual ~speedOfSound();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();

        //- Clear the result field from the objectRegistry
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const speedOfSound&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
