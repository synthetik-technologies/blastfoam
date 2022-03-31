/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    hydrostaticPressure

Description
    FunctionObject that outputs the hydrostatic component of the stress tensor
    field.

    Calculate and writes hydrostatic pressure field,
    where:
        hydrostaticPressure = -tr(sigma)/3.0;
        sigma is the stress tensor field.

    @verbatim
    hydrostaticPressure1
    {
        // Type of functionObject
        type hydrostaticPressure;

        // Where to load it from (if not already in solver)
        // Note: this is not required if the solver already loads this library.
        functionObjectLibs ("libsolidModelsBekaert.so");

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;
    );
    @endverbatim

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    hydrostaticPressure.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticPressure_H
#define hydrostaticPressure_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class hydrostaticPressure Declaration
\*---------------------------------------------------------------------------*/

class hydrostaticPressure
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        hydrostaticPressure
        (
            const hydrostaticPressure&
        );

        //- Disallow default bitwise assignment
        void operator=(const hydrostaticPressure&);


public:

    //- Runtime type information
    TypeName("hydrostaticPressure");


    // Constructors

        //- Construct from components
        hydrostaticPressure
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Destructor

        virtual ~hydrostaticPressure();


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
