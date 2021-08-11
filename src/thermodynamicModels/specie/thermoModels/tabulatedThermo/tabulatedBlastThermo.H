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
    Foam::tabulatedThermo

Description
    Thermodynamic model that uses a table based on density and internal energy
    to return temperature. Numerical derivatives are used to approximate
    additional quantities such as specific heats

SourceFiles
    tabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedThermoBlast_H
#define tabulatedThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTable1D.H"
#include "CpCvModel.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class tabulatedThermo;

template<class EquationOfState>
inline tabulatedThermo<EquationOfState> operator+
(
    const tabulatedThermo<EquationOfState>&,
    const tabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline tabulatedThermo<EquationOfState> operator*
(
    const scalar,
    const tabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline tabulatedThermo<EquationOfState> operator==
(
    const tabulatedThermo<EquationOfState>&,
    const tabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const tabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class tabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class tabulatedThermo
:
    public EquationOfState,
    public CpCvModel
{
// Private data

    //- Internal energy lookup table
    lookupTable1D eTable_;

    //- Temperature bounds
    scalar Tlow_;
    scalar Thigh_;


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

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


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


    // Member operators

        inline void operator+=(const tabulatedThermo&);


    // Friend operators

        friend tabulatedThermo operator+ <EquationOfState>
        (
            const tabulatedThermo&,
            const tabulatedThermo&
        );

        friend tabulatedThermo operator* <EquationOfState>
        (
            const scalar,
            const tabulatedThermo&
        );

        friend tabulatedThermo operator== <EquationOfState>
        (
            const tabulatedThermo&,
            const tabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedBlastThermoI.H"

#ifdef NoRepository
    #include "tabulatedBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
