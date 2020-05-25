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
    Foam::tabulated

Description
    Thermodynamic model that uses a table based on density and internal energy
    to return temperature. Numerical derivatives are used to approximate
    additional quantities such as specific heats

SourceFiles
    tabulated.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTable2D.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class tabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class tabulatedThermo
:
    public EquationOfState
{
// Private data

    //- Temperature lookup table
    lookupTable2D TTable_;


protected:

    // Protected data

        //- Convergence criteria for temperature
        scalar tolerance_;

        //- Maximum number of iterations
        label maxIter_;


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "tabulated<"
              + word(EquationOfState::typeName()) + ">";
        }

    // Constructors

        tabulatedThermo(const dictionary& dict);

    // Member Functions


        // Thermodynamic functions

            //- Return Mie Gruniesen coefficient
            scalar Gamma
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return speed of sound
            scalar speedOfSound
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            scalar Cv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            scalar Cp
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return specific heat ratio
            scalar CpByCv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Difference in specific heat ratios
            scalar CpMCv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Sensible internal energy [J/kg]
            scalar Es
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Absolute internal energy [J/kg]
            scalar Ea
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Sensible enthalpy [J/kg]
            scalar Hs
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const
            {
                return Cp(rho, e, T)*T + EquationOfState::H(rho, e, T);
            }

            //- Absolute enthalpy [J/kg]
            scalar Ha
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const
            {
                return Hs(rho, e, T);
            }

            //- Entropy [J/kg/K]
            scalar S
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return temperature
            scalar TRhoE
            (
                const scalar& T0,
                const scalar& rho,
                const scalar& e
            ) const;

            //- Initialize internal energy
            scalar initializeEnergy
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Initialize internal energy
            scalar initializeRho
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedThermoI.H"

#ifdef NoRepository
    #include "tabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
