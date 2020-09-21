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
    Foam::interfacialVelocityModels::single

Description
    One phase is used to determine the interfacial Velocity

SourceFiles
    singleInterfacialVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef singleInterfacialVelocity_H
#define singleInterfacialVelocity_H

#include "interfacialVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfacialVelocityModels
{

/*---------------------------------------------------------------------------*\
                          Class single Declaration
\*---------------------------------------------------------------------------*/

class single
:
    public interfacialVelocityModel
{
    //- Private data

        //- Name of phase
        const word phaseName_;

        //- const reference to phase model
        const phaseModel& phase_;

        //- Interfacial Velocity
        const volVectorField& U_;

        //- Interfacial volumeric flux
        const surfaceScalarField& phi_;


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

        //- Return the interfacial velocity
        virtual tmp<volVectorField> Ui() const;

        //- Return the interfacial volumetric flux
        virtual tmp<surfaceScalarField> phi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfacialVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
