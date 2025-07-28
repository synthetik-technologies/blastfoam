/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025
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
    Foam::functionObjects::stopProbes

Description
    Check field values at locations and stops simulation when the perscribed
    value is met

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    stopProbes.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stopProbes_H
#define functionObjects_stopProbes_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class stopProbes Declaration
\*---------------------------------------------------------------------------*/

class stopProbes
:
    public fvMeshFunctionObject
{
    // Private Member data

        //- Table of field name and values to stop at
        HashTable<scalar> fieldValues_;

        //- Sample locations
        pointField locations_;

        //- Stop control
        Time::stopAtControl stopControl_;

        //- Does the simulation need to be stopped
        bool stop_;

public:

    //- Runtime type information
    TypeName("stopProbes");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        stopProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        stopProbes(const stopProbes&) = delete;


    //- Destructor
    virtual ~stopProbes();




    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Check if the simulation should be stopped
        virtual bool execute();

        //- Write
        virtual bool write()
        {
            return true;
        }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
