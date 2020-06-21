/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Added writeTimeList functionObject
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
    Foam::functionObjects::writeTimeList

Description
    Overrides write interval to write at a designated list of times. The normal
    write interval specified in the control dict is still used.

    Example of function object specification:
    \verbatim
    max
    {
        type                writeTimeList;
        libs                ("libblastFunctionObjects.so");

        times (0.001 0.004 0.01);
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        time              | List of times to write at |yes
    \endtable


See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    writeTimeList.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeTimeList_H
#define functionObjects_writeTimeList_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class writeTimeList Declaration
\*---------------------------------------------------------------------------*/

class writeTimeList
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- List of write times
        scalarList writeTimes_;

        //- Index of next writeTime
        label index_;


public:

    //- Runtime type information
    TypeName("writeTimeList");


    // Constructors

        //- Construct from Time and dictionary
        writeTimeList
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeTimeList(const writeTimeList&) = delete;


    //- Destructor
    virtual ~writeTimeList();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Update index
        virtual bool execute();

        //- Return time to next write
        virtual scalar timeToNextWrite();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeTimeList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
