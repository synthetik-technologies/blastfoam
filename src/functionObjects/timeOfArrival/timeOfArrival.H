/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
20-06-2020 Synthetik Applied Technologies: | Time of arrival implementation
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
    Foam::functionObjects::timeOfArrival

Description
    Calculates and writes the time of arrival

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    timeOfArrival.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_timeOfArrival_H
#define functionObjects_timeOfArrival_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class timeOfArrival Declaration
\*---------------------------------------------------------------------------*/

class timeOfArrival
:
    public fvMeshFunctionObject
{
    // Private Member data

        //- Name of pressure field
        word pName_;

        //- Time of arrival field
        volScalarField timeOfArrival_;

        //- Max pressure field
        volScalarField& pMax_;

        //- Return the given field or create
        volScalarField& lookupOrCreate
        (
            const word& name,
            const dimensionSet& dims
        ) const;


public:

    //- Runtime type information
    TypeName("timeOfArrival");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        timeOfArrival
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        timeOfArrival(const timeOfArrival&) = delete;


    //- Destructor
    virtual ~timeOfArrival();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const timeOfArrival&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
