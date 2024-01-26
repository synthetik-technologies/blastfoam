/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate thermalImpulse
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
    Foam::functionObjects::thermalImpulse

Description
    Calculates the exposure due to temperature

    Example of function object specification:
    \verbatim
    thermalImpulse
    {
        type                thermalImpulse;
        libs                ("libblastFunctionObjects.so");

        writeControl    writeTime;
        restartOnRestart    false;

        TName           T;
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        restartOnRestart  | Restart the averaging on restart |no | no
        TName             | Name of temperature field | no       | T
    \endtable


See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    thermalImpulse.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_thermalImpulse_H
#define functionObjects_thermalImpulse_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class thermalImpulse Declaration
\*---------------------------------------------------------------------------*/

class thermalImpulse
:
    public fvMeshFunctionObject
{
public:
    enum class burnModel
    {
        CEM43,
        HENRIQUES
    };
    static const NamedEnum<burnModel, 2> burnModelNames_;

protected:

    // Protected data

        //- Restart the averaging process on output
        Switch restartOnRestart_;

        //- Name of temperature field
        word TName_;

        //- Name of qr (radiation heat flux) field
        word qrName_;

        //- Integral of temperature gradient
        // averaged snGrad
        volScalarField intGradT_;

        //- Integral of external energy flux (i.e. radiation heat flux)
        volScalarField intQExt_;

public:

    //- Runtime type information
    TypeName("thermalImpulse");


    // Constructors

        //- Construct from Time and dictionary
        thermalImpulse
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        thermalImpulse(const thermalImpulse&) = delete;


    //- Destructor
    virtual ~thermalImpulse();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermalImpulse&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
