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
    Foam::hConst

Description
    Thermodynamics model based with constant properties

SourceFiles
    hConstI.H
    hConst.C
\*---------------------------------------------------------------------------*/

#ifndef hConstBlast_H
#define hConstBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModel.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hConst;

template<class EquationOfState>
inline hConst<EquationOfState> operator+
(
    const hConst<EquationOfState>&,
    const hConst<EquationOfState>&
);

template<class EquationOfState>
inline hConst<EquationOfState> operator*
(
    const scalar,
    const hConst<EquationOfState>&
);

template<class EquationOfState>
inline hConst<EquationOfState> operator==
(
    const hConst<EquationOfState>&,
    const hConst<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hConst<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hConst Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hConst
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp_;

        //- Heat of formation [J/kg]
        scalar Hf_;


        //- Construct from components
        inline hConst
        (
            const EquationOfState& st,
            const scalar cv,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        hConst(const dictionary& dict);

        //- Construct as named copy
        inline hConst(const word&, const hConst&);

        //- Construct and return a clone
        inline autoPtr<hConst> clone() const;

        // Selector from dictionary
        inline static autoPtr<hConst> New(const dictionary& dict);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "hConst<"
              + word(EquationOfState::typeName()) + '>';
        }


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

        inline void operator+=(const hConst&);


    // Friend operators

        friend hConst operator+ <EquationOfState>
        (
            const hConst&,
            const hConst&
        );

        friend hConst operator* <EquationOfState>
        (
            const scalar,
            const hConst&
        );

        friend hConst operator== <EquationOfState>
        (
            const hConst&,
            const hConst&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hConstI.H"

#ifdef NoRepository
    #include "hConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
