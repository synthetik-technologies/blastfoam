/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::interfacialPressureModels::single

Description
    One phase is used to determine the interfacial pressure

SourceFiles
    singleInterfacialPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef singleInterfacialPressure_H
#define singleInterfacialPressure_H

#include "interfacialPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfacialPressureModels
{

/*---------------------------------------------------------------------------*\
                          Class single Declaration
\*---------------------------------------------------------------------------*/

class single
:
    public interfacialPressureModel
{
    //- Private data

        //- Name of phase
        const word phaseName_;

        //- const reference to phase model
        const phaseModel& phase_;

        //- Interfacial pressure
        const volScalarField& p_;


public:

    //- Runtime type information
    TypeName("single");


    // Constructors

        //- Construct from components
        single
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~single();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        virtual tmp<volScalarField> Pi() const;

        //- Return the interfacial pressure for a cell
        virtual scalar Pi(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfacialPressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
