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
    Foam::interfacialPressureModels::totalPressure

Description
    Includes thermodynamic and kinetic energy contributions

SourceFiles
    totalPressurePressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureInterfacialPressure_H
#define totalPressureInterfacialPressure_H

#include "interfacialPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfacialPressureModels
{

/*---------------------------------------------------------------------------*\
                          Class totalPressure Declaration
\*---------------------------------------------------------------------------*/

class totalPressure
:
    public interfacialPressureModel
{
    //- Reference to the interfacial velocity
    const volVectorField& U_;


public:

    //- Runtime type information
    TypeName("totalPressure");


    // Constructors

        //- Construct from components
        totalPressure
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~totalPressure();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        virtual tmp<volScalarField> PI() const;

        //- Return the interfacial pressure for a cell
        virtual scalar cellPI(const label celli) const;

        //- Return the derivative of the given phase w.r.t. volume fraction
        virtual scalar celldPIdAlpha
        (
            const label celli,
            const label phasei
        ) const;

        //- Return the derivative of the given phase w.r.t. internal energy
        virtual scalar celldPIde
        (
            const label celli,
            const label phasei
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfacialPressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
