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
    Foam::janaf9Thermo

Description
    janaf9Thermo tables based thermodynamics package templated
    into the equation of state.

    entries 0:n-2 are cp coefficients, n-1 is the enthaply offset, and n-2 is the entropy offset

SourceFiles
    janaf9ThermoBlastThermo.C
\*---------------------------------------------------------------------------*/

#ifndef janaf9ThermoBlast_H
#define janaf9ThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "List.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class janaf9Thermo;

template<class EquationOfState>
inline janaf9Thermo<EquationOfState> operator+
(
    const janaf9Thermo<EquationOfState>&,
    const janaf9Thermo<EquationOfState>&
);

template<class EquationOfState>
inline janaf9Thermo<EquationOfState> operator*
(
    const scalar,
    const janaf9Thermo<EquationOfState>&
);

template<class EquationOfState>
inline janaf9Thermo<EquationOfState> operator==
(
    const janaf9Thermo<EquationOfState>&,
    const janaf9Thermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const janaf9Thermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class janaf9Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class janaf9Thermo
:
    public EquationOfState
{
public:

    // Public data

        static const int nCoeffs_ = 9;
        typedef FixedList<scalar, nCoeffs_> coeffArray;



private:

    // Private Data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private Member Functions

        //- Return the coefficients corresponding to the given temperature
        const coeffArray& coeffs(const scalar T) const;



public:

    // Constructors

        //- Construct from components
        inline janaf9Thermo
        (
            const EquationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon,
            const coeffArray& highCpCoeffs,
            const coeffArray& lowCpCoeffs,
            const bool convertCoeffs = false
        );

        //- Construct from dictionary
        janaf9Thermo(const dictionary& dict);

        //- Construct as a named copy
        inline janaf9Thermo(const word&, const janaf9Thermo&);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "janaf9<"
              + word(EquationOfState::typeName()) + '>';
        }

        //- Return the energy type
        static bool enthalpy()
        {
            return true;
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T, const scalar Told) const;


    // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the common temperature
            inline scalar Tcommon() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& highCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& lowCpCoeffs() const;


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

        inline void operator+=(const janaf9Thermo&);


    // Friend operators

        friend janaf9Thermo operator+ <EquationOfState>
        (
            const janaf9Thermo&,
            const janaf9Thermo&
        );

        friend janaf9Thermo operator* <EquationOfState>
        (
            const scalar,
            const janaf9Thermo&
        );

        friend janaf9Thermo operator== <EquationOfState>
        (
            const janaf9Thermo&,
            const janaf9Thermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const janaf9Thermo&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "janaf9BlastThermoI.H"

#ifdef NoRepository
    #include "janaf9BlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
