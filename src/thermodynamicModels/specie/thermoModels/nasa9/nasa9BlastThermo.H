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
    Foam::nasa9Thermo

Description
    NASA 9-coefficient polynomial thermodynamics (McBride, Zehe and Gordon,
    NASA/TP-2002-211556) with any number of temperature ranges, templated
    into the equation of state.

    In each range
        Cp/R = a0/T^2 + a1/T + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
        H/R  = -a0/T + a1 ln(T) + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4
             + a6 T^5/5 + b1
        S/R  = -a0/(2T^2) - a1/T + a2 ln(T) + a3 T + a4 T^2/2 + a5 T^3/3
             + a6 T^4/4 + b2
    with the coefficients (a0 .. a6 b1 b2) per range. Outside the overall
    range Tswitch[0] - Tswitch[n] the properties are extrapolated with
    constant Cp (H linear in T, S extended with Cp ln(T/Tlimit)).

    Mixtures are formed by mass-weighting the coefficients, so all species
    of a mixture must use the same temperature ranges.

    \verbatim
        thermodynamics
        {
            Tswitch     (200 1000 6000 20000);
            coeffs
            (
                (a0 a1 a2 a3 a4 a5 a6 b1 b2)    // 200 - 1000 K
                (...)                           // 1000 - 6000 K
                (...)                           // 6000 - 20000 K
            );
        }
    \endverbatim

SourceFiles
    nasa9BlastThermoI.H
    nasa9BlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef nasa9ThermoBlast_H
#define nasa9ThermoBlast_H

#include "scalar.H"
#include "FixedList.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class nasa9Thermo;

template<class EquationOfState>
inline nasa9Thermo<EquationOfState> operator+
(
    const nasa9Thermo<EquationOfState>&,
    const nasa9Thermo<EquationOfState>&
);

template<class EquationOfState>
inline nasa9Thermo<EquationOfState> operator*
(
    const scalar,
    const nasa9Thermo<EquationOfState>&
);

template<class EquationOfState>
inline nasa9Thermo<EquationOfState> operator==
(
    const nasa9Thermo<EquationOfState>&,
    const nasa9Thermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const nasa9Thermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                         Class nasa9Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class nasa9Thermo
:
    public EquationOfState
{
public:

    // Public data

        static const int nCoeffs_ = 9;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private Data

        //- Temperature range limits (n + 1 values for n ranges)
        List<scalar> Tswitch_;

        //- Coefficients per range, on a mass basis (multiplied by R)
        List<coeffArray> coeffs_;


    // Private Member Functions

        //- Index of the range containing T (T within the overall range)
        inline label range(const scalar T) const;

        //- Return the temperature clipped to the overall range
        inline scalar Tfit(const scalar T) const;

        //- Polynomial Cp/H/S (standard state) at a temperature in the range
        inline scalar CpFit(const scalar Tf) const;
        inline scalar HaFit(const scalar Tf) const;
        inline scalar SFit(const scalar Tf) const;

        //- Standard-state entropy, extrapolated with constant Cp outside
        //  the range
        inline scalar Sstd(const scalar T) const;

        //- Check that two thermos have the same temperature ranges
        static void checkRanges(const nasa9Thermo&, const nasa9Thermo&);


public:

    // Constructors

        //- Construct from components
        inline nasa9Thermo
        (
            const EquationOfState& st,
            const List<scalar>& Tswitch,
            const List<coeffArray>& coeffs,
            const bool convertCoeffs = false
        );

        //- Construct from dictionary
        nasa9Thermo(const dictionary& dict);

        //- Construct as a named copy
        inline nasa9Thermo(const word&, const nasa9Thermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "nasa9<"
              + word(EquationOfState::typeName()) + '>';
        }

        //- Return the energy type
        static bool enthalpy()
        {
            return true;
        }

        //- Limit the temperature to Tlow to 10*Thigh (above Thigh the
        //  properties are extrapolated with constant Cp)
        inline scalar limit(const scalar T, const scalar Told) const;


        // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the range limits
            inline const List<scalar>& Tswitch() const;

            //- Return const access to the coefficients (mass basis)
            inline const List<coeffArray>& coeffs() const;


        // Fundamental properties

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

            //- Absolute enthalpy [J/kg]
            inline scalar Ha
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Heat capacity at constant pressure of mixture [J/kg/K]
            inline scalar Cpt(const scalar T) const;

            //- Sensible thermal enthalpy [J/kg]
            inline scalar Hst(const scalar T) const;

            //- Absolute thermal enthalpy [J/kg]
            inline scalar Hat(const scalar T) const;

            #include "HtoEFuncs.H"

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


    // Member Operators

        inline void operator+=(const nasa9Thermo&);


    // Friend operators

        friend nasa9Thermo operator+ <EquationOfState>
        (
            const nasa9Thermo&,
            const nasa9Thermo&
        );

        friend nasa9Thermo operator* <EquationOfState>
        (
            const scalar,
            const nasa9Thermo&
        );

        friend nasa9Thermo operator== <EquationOfState>
        (
            const nasa9Thermo&,
            const nasa9Thermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const nasa9Thermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nasa9BlastThermoI.H"

#ifdef NoRepository
    #include "nasa9BlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
