/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate overpressure
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
    Foam::functionObjects::overpressure

Description
    Calculates and writes the Mach number as a volScalarField.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    overpressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_overpressure_H
#define functionObjects_overpressure_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class overpressure Declaration
\*---------------------------------------------------------------------------*/

class overpressure
:
    public fvMeshFunctionObject
{
    // Private Member data

        //- Name of pressure
        word pName_;

        //- Name of Result field
        word resultName_;

        //- Reference pressure
        dimensionedScalar pRef_;

        //- Is the field stored in the database
        Switch store_;


public:

    //- Runtime type information
    TypeName("overpressure");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        overpressure
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        overpressure(const overpressure&) = delete;


    //- Destructor
    virtual ~overpressure();


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
        void operator=(const overpressure&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
