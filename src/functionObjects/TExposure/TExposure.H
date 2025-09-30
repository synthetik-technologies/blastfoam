/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate TExposure
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
    Foam::functionObjects::TExposure

Description
    Calculates the exposure due to temperature

    Example of function object specification:
    \verbatim
    TExposure
    {
        type                TExposure;
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
    TExposure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_TExposure_H
#define functionObjects_TExposure_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class TExposure Declaration
\*---------------------------------------------------------------------------*/

class TExposure
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

        //- Exposure field
        volScalarField exposure_;

        //- Burn model
        burnModel burnModel_;

        //- Pre exponential factor
        Pair<scalar> P_;

        //- Activation energy ratio
        Pair<scalar> EbyR_;

        //- Switching temperature
        scalar Tswitch_;

        //- Minimum temperature
        scalar Tmin_;

        //- Maximum temperature
        scalar Tmax_;

public:

    //- Runtime type information
    TypeName("TExposure");


    // Constructors

        //- Construct from Time and dictionary
        TExposure
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        TExposure(const TExposure&) = delete;


    //- Destructor
    virtual ~TExposure();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return {TName_};
        }

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TExposure&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
