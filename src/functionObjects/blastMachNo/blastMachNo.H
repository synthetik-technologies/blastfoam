/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate mach number with blastFoam thermo
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
    Foam::functionObjects::MachNo

Description
    Calculates and writes the Mach number as a volScalarField.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    blastMachNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blastMachNo_H
#define functionObjects_blastMachNo_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class MachNo Declaration
\*---------------------------------------------------------------------------*/

template<class CompressibleSystem>
class blastMachNo
:
    public fvMeshFunctionObject
{
    // Private Member data

        static const word systemType_;

        //- Name of phase
        word phaseName_;

        //- Name of looked up system
        word systemName_;

        //- Name of Result field
        word resultName_;

        //- Name of velocity field
        word UName_;


public:


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        blastMachNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        blastMachNo(const blastMachNo&) = delete;


    //- Destructor
    virtual ~blastMachNo();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();

        //- Clear the result field from the objectRegistry
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastMachNo&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMachNo(systemType, name)                                        \
template<class systemType>                                                  \
const word blastMachNo<systemType>::systemType_ = #systemType;              \
                                                                            \
typedef blastMachNo<systemType> systemType##MachNo;                         \
                                                                            \
addNamedToRunTimeSelectionTable                                             \
(                                                                           \
    functionObject,                                                         \
    systemType##MachNo,                                                     \
    dictionary,                                                             \
    name                                                                    \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blastMachNo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
