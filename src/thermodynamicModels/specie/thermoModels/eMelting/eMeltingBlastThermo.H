/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::eMeltingThermo

Description
    Internal-energy based thermodynamics of a material that melts: constant
    heat capacity of the solid (Cv) below and of the liquid (CvLiquid) above
    the melting temperature Tm, and the latent heat of melting Lm released
    over a narrow range dTm centred on Tm (the heat capacity there is
    Cv + Lm/dTm, so that it stays finite):

        T1 = Tm - dTm/2,  T2 = Tm + dTm/2
        Es(T) = Esref + Cv (min(T, T1) - Tref)
              + (Cv + Lm/dTm) min(max(T - T1, 0), dTm)
              + CvLiquid max(T - T2, 0)
              + E_eos(rho)

    Tref must be below T1. The energy is piecewise linear in T and is
    inverted exactly (TEsEstimate), since Newton's method on the kinked
    curve can cycle between the solid and liquid branches.

    \verbatim
    thermodynamics
    {
        type        eMelting;
        Cv          1040;       // Solid, mean over 300 - 933 K [J/kg/K]
        CvLiquid    1177;       // Liquid [J/kg/K]
        Tm          933.47;     // Melting temperature [K]
        Lm          3.97e5;     // Latent heat of melting [J/kg]
        dTm         2;          // Melting range [K] (optional, default 2)
        Hf          0;
        Tref        298.15;     // Optional, default Tstd
    }
    \endverbatim

SourceFiles
    eMeltingBlastThermoI.H
    eMeltingBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eMeltingBlastThermo_H
#define eMeltingBlastThermo_H

#include "thermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eMeltingThermo;

template<class EquationOfState>
inline eMeltingThermo<EquationOfState> operator+
(
    const eMeltingThermo<EquationOfState>&,
    const eMeltingThermo<EquationOfState>&
);

template<class EquationOfState>
inline eMeltingThermo<EquationOfState> operator*
(
    const scalar,
    const eMeltingThermo<EquationOfState>&
);

template<class EquationOfState>
inline eMeltingThermo<EquationOfState> operator==
(
    const eMeltingThermo<EquationOfState>&,
    const eMeltingThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eMeltingThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                       Class eMeltingThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eMeltingThermo
:
    public EquationOfState
{
    // Private data

        //- Heat capacity of the solid [J/kg/K]
        scalar Cv_;

        //- Heat capacity of the liquid [J/kg/K]
        scalar CvLiquid_;

        //- Melting temperature [K]
        scalar Tm_;

        //- Latent heat of melting [J/kg]
        scalar Lm_;

        //- Temperature range of melting [K]
        scalar dTm_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Reference temperature [K]
        scalar Tref_;

        //- Sensible energy at Tref [J/kg]
        scalar Esref_;


    // Private Member Functions

        //- Start of melting [K]
        inline scalar T1() const;

        //- End of melting [K]
        inline scalar T2() const;

        //- Thermal energy relative to Tref [J/kg]
        inline scalar dEt(const scalar T) const;

        //- Construct from components
        inline eMeltingThermo
        (
            const EquationOfState& st,
            const scalar Cv,
            const scalar CvLiquid,
            const scalar Tm,
            const scalar Lm,
            const scalar dTm,
            const scalar Hf,
            const scalar Tref,
            const scalar Esref
        );


public:

    // Constructors

        //- Construct from dictionary
        eMeltingThermo(const dictionary& dict);

        //- Construct as named copy
        inline eMeltingThermo(const word&, const eMeltingThermo&);

        //- Construct and return a clone
        inline autoPtr<eMeltingThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<eMeltingThermo> New(const dictionary& dict);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "eMelting<"
              + word(EquationOfState::typeName()) + '>';
        }

        //- Return the energy type
        static bool enthalpy()
        {
            return false;
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T, const scalar) const;


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

        //- Absolute internal energy [J/kg]
        inline scalar Ea
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Temperature from the sensible internal energy (exact inverse
        //  of Es; T0 is used for the equation of state energy)
        inline scalar TEsEstimate
        (
            const scalar rho,
            const scalar es,
            const scalar T0
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        inline scalar Cvt(const scalar T) const;

        //- Sensible thermal internal energy [J/kg]
        inline scalar Est(const scalar T) const;

        //- Absolute thermal internal energy [J/kg]
        inline scalar Eat(const scalar T) const;

        #include "EtoHFuncs.H"

        //- Heat of formation [J/kg]
        inline scalar Hf() const;

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


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const eMeltingThermo&);


    // Friend operators

        friend eMeltingThermo operator+ <EquationOfState>
        (
            const eMeltingThermo&,
            const eMeltingThermo&
        );

        friend eMeltingThermo operator* <EquationOfState>
        (
            const scalar,
            const eMeltingThermo&
        );

        friend eMeltingThermo operator== <EquationOfState>
        (
            const eMeltingThermo&,
            const eMeltingThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const eMeltingThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eMeltingBlastThermoI.H"

#ifdef NoRepository
    #include "eMeltingBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
