/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate impulse
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
    Foam::functionObjects::impulse

Description
    Calculates the impulse using the difference from a reference pressure

    Example of function object specification:
    \verbatim
    impulse
    {
        type                impulse;
        libs                ("libblastFunctionObjects.so");

        writeControl    writeTime;
        restartOnRestart    false;

        pName           p;
        pRef            101298;
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        type              | type name: fieldAverage   | yes      |
        restartOnRestart  | Restart the averaging on restart |no | no
        pName             | Name of pressure field    | no       | p
        pRef              | Reference pressure        | yes      |
    \endtable


See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    impulse.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_impulse_H
#define functionObjects_impulse_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class impulse Declaration
\*---------------------------------------------------------------------------*/

class impulse
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Restart the averaging process on output
        Switch restartOnRestart_;

        //- Name of pressure field
        const volScalarField& p_;

        //- Reference pressure
        dimensionedScalar pRef_;

        //- Saved impulse field
        volScalarField& impulse_;

        //- Create impulse field and add to object registry
        volScalarField& createImpulseField(const word& name);

public:

    //- Runtime type information
    TypeName("impulse");


    // Constructors

        //- Construct from Time and dictionary
        impulse
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        impulse(const impulse&) = delete;


    //- Destructor
    virtual ~impulse();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const impulse&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
