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
    Foam::janaf

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    janafThermo.C
\*---------------------------------------------------------------------------*/

#ifndef janafBlast_H
#define janafBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "FixedList.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class janaf;

template<class EquationOfState>
inline janaf<EquationOfState> operator+
(
    const janaf<EquationOfState>&,
    const janaf<EquationOfState>&
);

template<class EquationOfState>
inline janaf<EquationOfState> operator*
(
    const scalar,
    const janaf<EquationOfState>&
);

template<class EquationOfState>
inline janaf<EquationOfState> operator==
(
    const janaf<EquationOfState>&,
    const janaf<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const janaf<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class janaf Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class janaf
:
    public EquationOfState
{
public:

    // Public data

        typedef List<scalar> coeffArray;


private:

    // Private Data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private Member Functions

        //- Return the coefficients corresponding to the given temperature
        const coeffArray& coeffs(const scalar& T) const;



public:

    // Constructors

        //- Construct from components
        inline janaf
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
        janaf(const dictionary& dict);

        //- Construct as a named copy
        inline janaf(const word&, const janaf&);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "janaf<"
              + word(EquationOfState::typeName()) + '>';
        }


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
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Sensible enthalpy [J/kg]
        inline scalar Hs
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Absolute enthalpy [J/kg]
        inline scalar Ha
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        #include "HtoEFuncs.H"

        //- Heat of formation [J/kg]
        inline scalar Hf() const;

        //- Entropy [J/kg/K]
        inline scalar S
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;


    // Member operators

        inline void operator+=(const janaf&);


    // Friend operators

        friend janaf operator+ <EquationOfState>
        (
            const janaf&,
            const janaf&
        );

        friend janaf operator* <EquationOfState>
        (
            const scalar,
            const janaf&
        );

        friend janaf operator== <EquationOfState>
        (
            const janaf&,
            const janaf&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "janafI.H"

#ifdef NoRepository
    #include "janaf.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
