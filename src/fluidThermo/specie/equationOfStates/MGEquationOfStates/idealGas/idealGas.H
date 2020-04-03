/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::idealGas

Description
    Equation of state for an ideal gas using the Mie Grunsien form

SourceFiles
    idealGas.C

\*---------------------------------------------------------------------------*/

#ifndef idealGasBlast_H
#define idealGasBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class idealGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class idealGas
:
    public Specie
{
// Private data

    //- Specific heat ratio
    scalar gamma_;


public:

    // Constructors

        //- Single phase constructor
        idealGas(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "idealGas<"
              + word(Specie::typeName()) + '>';
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return false;
        }

        //- Return deviation from ideal gas pressure
        scalar Pi
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Return non-ideal speed of sound
        scalar delta
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e
        ) const;

        //- Gamma function
        scalar Gamma
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpde
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar& p,
            const scalar& rho,
            const scalar& T
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "idealGasI.H"

#ifdef NoRepository
    #include "idealGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
