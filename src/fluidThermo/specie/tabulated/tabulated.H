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
    Equation of state that uses a table to return equation of state variables

SourceFiles
    tabulated.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTable.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tabulated Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class tabulated
:
    public Specie
{
// Private data

    //- Pressure lookup table
    lookupTable pTable_;

    //- Temperature lookup table
    lookupTable TTable_;


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "tabulated<tabulated<"
              + word(Specie::typeName()) + ">>";
        }

    // Constructors

        tabulated(const dictionary& dict);

    // Member Functions

        //- Is temperature used
        bool temperatureBased() const
        {
            return false;
        }

        // Equation of state functions

            //- Return pressure
            scalar p(const scalar& rho, const scalar& e, const scalar& T) const;

            //- Return derivative of pressure w.r.t. specific volume
            scalar dpdv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return derivative of pressure w.r.t. Temperature
            scalar dpdT
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return derivative of pressure w.r.t. internal energy
            scalar dpde
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

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

            //- Heat capacity difference [J/kg/K]
            scalar CpMCv
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
            ) const;

            //- Absolute enthalpy [J/kg]
            scalar Ha
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

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

#include "tabulatedI.H"

#ifdef NoRepository
    #include "tabulated.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
