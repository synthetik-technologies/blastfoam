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
    Foam::hTabulatedThermo

Description
    Thermodynamic model that uses a table based on density and internal energy
    to return temperature. Numerical derivatives are used to approximate
    additional quantities such as specific heats

SourceFiles
    hTabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hTabulatedThermoBlast_H
#define hTabulatedThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTables1D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hTabulatedThermo;

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class hTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hTabulatedThermo
:
    public EquationOfState
{
// Private data

    //- Enthalpy lookup table
    scalarLookupTable1D hTable_;

    //- Temperature bounds
    scalar Tlow_;
    scalar Thigh_;


    //- Construct from components
    inline hTabulatedThermo
    (
        const EquationOfState& st,
        const scalarLookupTable1D& hTable
    );


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "hTabulated<"
              + word(EquationOfState::typeName()) + ">";
        }

    // Constructors

        hTabulatedThermo(const dictionary& dict);

        //- Return the energy type
        static bool enthalpy()
        {
            return true;
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T, const scalar Told) const;


    // Member Functions

        //- Heat capacity at constant pressure of mixture [J/kg/K]
        inline scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Sensible enthalpy [J/kg]
        inline scalar Hs
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Sensible thermal enthalpy [J/kg]
        inline scalar Hts
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Absolute enthalpy [J/kg]
        inline scalar Ha
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Absolute thermal enthalpy [J/kg]
        inline scalar Hta
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        #include "HtoEFuncs.H"

        //- Heat of formation [J/kg]
        inline scalar Hf() const;

        //- Flame temperature [K]
        inline scalar flameT() const;

        //- Entropy [J/kg/K]
        inline scalar S
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Derivative of Gibbs free energy w.r.t. temperature
        inline scalar dGdT(const scalar p, const scalar T) const;

        //- Temperature derivative of heat capacity at constant pressure
        inline scalar dCpdT(const scalar p, const scalar T) const;

        //- Temperature derivative of heat capacity at constant pressure
        inline scalar dCpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Temperature derivative of heat capacity at constant volume
        inline scalar dCvdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTabulatedBlastThermoI.H"

#ifdef NoRepository
    #include "hTabulatedBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
