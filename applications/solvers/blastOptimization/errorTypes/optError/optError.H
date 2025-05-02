/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies: | Calculate optError
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
    Foam::functionObjects::optError

Description
    Calculates and writes the Mach number as a volScalarField.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    optError.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_optError_H
#define functionObjects_optError_H

#include "fvMeshFunctionObject.H"
#include "errorType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class optError Declaration
\*---------------------------------------------------------------------------*/

class optError
:
    public fvMeshFunctionObject
{
    // Private Member data

        //- Output file name
        fileName output_;

        //- Error criteria
        PtrList<errorType> errors_;


public:

    //- Runtime type information
    TypeName("optError");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        optError
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        optError(const optError&) = delete;


    //- Destructor
    virtual ~optError();


    // Member Functions

        //- Reuired Fields
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the error to the given file
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const optError&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
