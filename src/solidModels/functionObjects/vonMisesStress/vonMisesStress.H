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
    Foam::functionObjects::vonMisesStress

Description
    FunctionObject that outputs the vonMises stress

    @verbatim
    vonMisesStress1
    {
        // Type of functionObject
        type vonMisesStress;

        // Where to load it from (if not already in solver)
        // Note: this is not required if the solver already loads this library.
        functionObjectLibs ("libsolidModels.so");

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;
    );
    @endverbatim

SourceFiles
    vonMisesStress.C

\*---------------------------------------------------------------------------*/

#ifndef vonMisesStress_H
#define vonMisesStress_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class vonMisesStress Declaration
\*---------------------------------------------------------------------------*/

class vonMisesStress
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("vonMises");


    // Constructors

        //- Construct from components
        vonMisesStress
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return {"sigma"};
        }

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
