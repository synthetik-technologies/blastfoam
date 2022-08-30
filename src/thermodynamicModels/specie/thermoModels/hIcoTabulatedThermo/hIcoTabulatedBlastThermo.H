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
    Foam::hIcoTabulatedThermo

Description
    Thermodynamic model that uses a table based on internal energy
    to return temperature. Numerical derivatives are used to approximate
    additional quantities such as specific heats

SourceFiles
    hIcoTabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hIcoTabulatedThermoBlast_H
#define hIcoTabulatedThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTables1D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hIcoTabulatedThermo;

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hIcoTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class hIcoTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hIcoTabulatedThermo
:
    public EquationOfState
{
// Private data

    //- Enthalpy lookup table
    scalarLookupTable1D hTable_;

    //- Temperature bounds
    scalar Tlow_;
    scalar Thigh_;

    //- Heat of formation
    scalar Hf_;


    //- Construct from components
    inline hIcoTabulatedThermo
    (
        const EquationOfState& st,
        const scalarLookupTable1D& hTable,
        const scalar Hf
    );


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "hIcoTabulated<"
              + word(EquationOfState::typeName()) + ">";
        }

    // Constructors

        hIcoTabulatedThermo(const dictionary& dict);

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
        inline scalar S(const scalar p, const scalar T) const;

        //- Gibbs free energy of the mixture in the standard state [J/kg]
        inline scalar Gstd(const scalar T) const;

        //- Temperature derivative of heat capacity at constant pressure
        inline scalar dCpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hIcoTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hIcoTabulatedBlastThermoI.H"

#ifdef NoRepository
    #include "hIcoTabulatedBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
