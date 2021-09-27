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
    Foam::hConstThermo

Description
    Enthalpy based thermodynamics package using a constant heat capacity
    at constant pressure:

    \verbatim
        h = Cp*(T - Tref) + Hsref
    \endverbatim

Usage
    \table
        Property     | Description
        Cp           | Constant Heat capacity at constant pressure [J/kmol/K]
        Tref         | Reference temperature [K] (defaults to 0)
        Hsref        | Reference sensible enthalpy [J/kg] (defaults to 0)
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of hConstThermo for air:
    \verbatim
    thermodynamics
    {
        Cp          1007;
        Hf          0;
    }
    \endverbatim

SourceFiles
    hConstThermoI.H
    hConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hConstThermoBlast_H
#define hConstThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModel.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hConstThermo;

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator+
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator*
(
    const scalar,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline hConstThermo<EquationOfState> operator==
(
    const hConstThermo<EquationOfState>&,
    const hConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hConstThermo
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant pressure [J/kmol/K]
        scalar Cp_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Reference temperature around which to linearise [K]
        scalar Tref_;

        //- Reference sensible enthalpy around which to linearise [J/kg]
        scalar Hsref_;

        //- Flame temperature [K]
        scalar flameT_;


        //- Construct from components
        inline hConstThermo
        (
            const EquationOfState& st,
            const scalar Cp,
            const scalar Hf,
            const scalar Tref,
            const scalar Hsref,
            const scalar flameT
        );


public:

    // Constructors

        //- Construct from dictionary
        hConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline hConstThermo(const word&, const hConstThermo&);

        //- Construct and return a clone
        inline autoPtr<hConstThermo> clone() const;

        // Selector from dictionary
        inline static autoPtr<hConstThermo> New(const dictionary& dict);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "hConst<"
              + word(EquationOfState::typeName()) + '>';
        }

        //- Return the energy type
        static bool enthalpy()
        {
            return true;
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T, const scalar) const;


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


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hConstThermo&);


    // Friend operators

        friend hConstThermo operator+ <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );

        friend hConstThermo operator* <EquationOfState>
        (
            const scalar,
            const hConstThermo&
        );

        friend hConstThermo operator== <EquationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include"hConstBlastThermoI.H"

#ifdef NoRepository
    #include "hConstBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
