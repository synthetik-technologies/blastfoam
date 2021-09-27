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
    Foam::ePowerBlastThermo

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

    Example specification of ePowerThermo:
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
    ePowerBlastThermoI.H
    ePowerBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ePowerThermoBlast_H
#define ePowerThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class ePowerThermo;

template<class EquationOfState>
inline ePowerThermo<EquationOfState> operator+
(
    const ePowerThermo<EquationOfState>&,
    const ePowerThermo<EquationOfState>&
);

template<class EquationOfState>
inline ePowerThermo<EquationOfState> operator*
(
    const scalar,
    const ePowerThermo<EquationOfState>&
);

template<class EquationOfState>
inline ePowerThermo<EquationOfState> operator==
(
    const ePowerThermo<EquationOfState>&,
    const ePowerThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const ePowerThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class ePowerThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class ePowerThermo
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
        inline ePowerThermo
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
                "ePower<"
              + word(EquationOfState::typeName()) + ">";
        }

    // Constructors

        //- Construct from dictionary
        ePowerThermo(const dictionary& dict);

        //- Construct as a named copy
        inline ePowerThermo
        (
            const word&,
            const ePowerThermo&
        );

        //- Construct and return a clone
        inline autoPtr<ePowerThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<ePowerThermo> New(const dictionary& dict);

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


    // Member operators

        inline void operator+=(const ePowerThermo&);


    // Friend operators

        friend ePowerThermo operator+ <EquationOfState>
        (
            const ePowerThermo&,
            const ePowerThermo&
        );

        friend ePowerThermo operator* <EquationOfState>
        (
            const scalar,
            const ePowerThermo&
        );

        friend ePowerThermo operator== <EquationOfState>
        (
            const ePowerThermo&,
            const ePowerThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const ePowerThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ePowerBlastThermoI.H"

#ifdef NoRepository
    #include "ePowerBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
