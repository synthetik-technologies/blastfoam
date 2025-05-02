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
    stressTriaxiality

Description
    FunctionObject that outputs the stress triaxiality field:
    where: stressTriaxiality = hydrostaticStress/equivalentStress;

    @verbatim
    stressTriaxiality1
    {
        // Type of functionObject
        type stressTriaxiality;

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
    stressTriaxiality.C

\*---------------------------------------------------------------------------*/

#ifndef stressTriaxiality_H
#define stressTriaxiality_H

#include "fvMeshFunctionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class stressTriaxiality Declaration
\*---------------------------------------------------------------------------*/

class stressTriaxiality
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("stressTriaxiality");


    // Constructors

        //- Construct from components
        stressTriaxiality
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
            return {"sigma"};
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
