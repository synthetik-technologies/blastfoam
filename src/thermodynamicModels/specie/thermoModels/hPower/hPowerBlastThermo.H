/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2021-09-27 Jeff Heylmun:    Modified class for a density based thermodynamic
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
    Foam::hPowerBlastThermo

Description
    Internal energy based thermodynamics package using a power function of
    temperature for the constant heat capacity at constant volume which is
    particularly suitable for solids at low temperatures:

    \verbatim
        Cv = c0*pow(T/Tref, n0)
    \endverbatim

Usage
    \table
        Property     | Description
        c0           | Reference heat capacity at constant volume [J/kmol/K]
        n0           | Exponent of the power function
        Tref         | Reference temperature [K]
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of hPowerThermo:
    \verbatim
    thermodynamics
    {
        c0          230;
        Tref        470;
        n0          3;
        Hf          0;
    }
    \endverbatim

SourceFiles
    hPowerBlastThermoI.H
    hPowerBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hPowerThermoBlast_H
#define hPowerThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hPowerThermo;

template<class EquationOfState>
inline hPowerThermo<EquationOfState> operator+
(
    const hPowerThermo<EquationOfState>&,
    const hPowerThermo<EquationOfState>&
);

template<class EquationOfState>
inline hPowerThermo<EquationOfState> operator*
(
    const scalar,
    const hPowerThermo<EquationOfState>&
);

template<class EquationOfState>
inline hPowerThermo<EquationOfState> operator==
(
    const hPowerThermo<EquationOfState>&,
    const hPowerThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hPowerThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class hPowerThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hPowerThermo
:
    public EquationOfState
{
    // Private Data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline hPowerThermo
        (
            const EquationOfState& st,
            const scalar c0,
            const scalar n0,
            const scalar Tref,
            const scalar Hf
        );


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "hPower<"
              + word(EquationOfState::typeName()) + ">";
        }

    // Constructors

        //- Construct from dictionary
        hPowerThermo(const dictionary& dict);

        //- Construct as a named copy
        inline hPowerThermo
        (
            const word&,
            const hPowerThermo&
        );

        //- Construct and return a clone
        inline autoPtr<hPowerThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hPowerThermo> New(const dictionary& dict);

        //- Return the energy type
        static bool enthalpy()
        {
            return true;
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T, const scalar Told) const;


    // Member Functions

        //- Heat capacity at constant pressure [J/kg/K]
        inline scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Absolute enthalpy [J/kg]
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


    // Member operators

        inline void operator+=(const hPowerThermo&);


    // Friend operators

        friend hPowerThermo operator+ <EquationOfState>
        (
            const hPowerThermo&,
            const hPowerThermo&
        );

        friend hPowerThermo operator* <EquationOfState>
        (
            const scalar,
            const hPowerThermo&
        );

        friend hPowerThermo operator== <EquationOfState>
        (
            const hPowerThermo&,
            const hPowerThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hPowerThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hPowerBlastThermoI.H"

#ifdef NoRepository
    #include "hPowerBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
