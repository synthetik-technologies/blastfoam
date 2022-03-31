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
    principalStresses

Description
    FunctionObject that outputs the principal stress fields.

    Calculate and write principal stress fields.

    Three vector fields are created: sigmaMax, sigmaMid, sigmaMin.
    sigmaMax is the most positive/tensile principal stress multiplied by the
    corresponding principal direction;
    sigmaMid is the middle principal stress multiplied by the corresponding
    principal direction;
    sigmaMin is the most negative/compressive principal stress multiplied by the
    corresponding principal direction.

    @verbatim
    principalStresses1
    {
        // Type of functionObject
        type principalStresses;

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;
    );
    @endverbatim

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    principalStresses.C

\*---------------------------------------------------------------------------*/

#ifndef principalStresses_H
#define principalStresses_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class principalStresses Declaration
\*---------------------------------------------------------------------------*/

class principalStresses
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& runTime_;

        //- Reference to mesh region
        const fvMesh& mesh_;


    // Private Member Functions

        //- Calculate and set principal stress vectors
        void calculateEigenValues
        (
            const symmTensor& sigma,
            vector& sigmaMax,
            vector& sigmaMid,
            vector& sigmaMin
        );

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        principalStresses
        (
            const principalStresses&
        );

        //- Disallow default bitwise assignment
        void operator=(const principalStresses&);


public:

    //- Runtime type information
    TypeName("principalStresses");


    // Constructors

        //- Construct from components
        principalStresses
        (
            const word& name,
            const Time&,
            const dictionary&
        );


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
