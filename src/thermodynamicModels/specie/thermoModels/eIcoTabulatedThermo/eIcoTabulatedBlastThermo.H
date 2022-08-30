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
    Foam::eIcoTabulatedThermo

Description
    Thermodynamic model that uses a table based on density and internal energy
    to return temperature. Numerical derivatives are used to approximate
    additional quantities such as specific heats

SourceFiles
    eIcoTabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eIcoTabulatedThermoBlast_H
#define eIcoTabulatedThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTables1D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eIcoTabulatedThermo;

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eIcoTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class eIcoTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eIcoTabulatedThermo
:
    public EquationOfState
{
// Private data

    //- Internal energy lookup table
    scalarLookupTable1D eTable_;

    //- Temperature bounds
    scalar Tlow_;
    scalar Thigh_;

    //- Heat of formation
    scalar Hf_;

    //- Construct from components
    inline eIcoTabulatedThermo
    (
        const EquationOfState& st,
        const scalarLookupTable1D& eTable,
        const scalar Hf
    );


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "eIcoTabulated<"
              + word(EquationOfState::typeName()) + ">";
        }

    // Constructors

        eIcoTabulatedThermo(const dictionary& dict);

        //- Return the energy type
        static bool enthalpy()
        {
            return false;
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T, const scalar Told) const;


    // Member Functions

        //- Heat capacity at constant volume of mixture [J/kg/K]
        inline scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Sensible internal energy [J/kg]
        inline scalar Es
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Sensible thermal internal energy [J/kg]
        inline scalar Ets
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Absolute internal energy [J/kg]
        inline scalar Ea
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Absolute thermal internal energy [J/kg]
        inline scalar Eta
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        #include "EtoHFuncs.H"

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
            const eIcoTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eIcoTabulatedBlastThermoI.H"

#ifdef NoRepository
    #include "eIcoTabulatedBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
