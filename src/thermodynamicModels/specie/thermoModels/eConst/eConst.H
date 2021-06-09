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
    Foam::eConst

Description
    Thermodynamics model based with constant properties

SourceFiles
    eConstI.H
    eConst.C
\*---------------------------------------------------------------------------*/

#ifndef eConstBlast_H
#define eConstBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModel.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eConst;

template<class EquationOfState>
inline eConst<EquationOfState> operator+
(
    const eConst<EquationOfState>&,
    const eConst<EquationOfState>&
);

template<class EquationOfState>
inline eConst<EquationOfState> operator*
(
    const scalar,
    const eConst<EquationOfState>&
);

template<class EquationOfState>
inline eConst<EquationOfState> operator==
(
    const eConst<EquationOfState>&,
    const eConst<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eConst<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class eConst Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eConst
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Flame temperature [K]
        scalar flameT_;


        //- Construct from components
        inline eConst
        (
            const EquationOfState& st,
            const scalar cv,
            const scalar hf,
            const scalar ft
        );


public:

    // Constructors

        //- Construct from dictionary
        eConst(const dictionary& dict);

        //- Construct as named copy
        inline eConst(const word&, const eConst&);

        //- Construct and return a clone
        inline autoPtr<eConst> clone() const;

        // Selector from dictionary
        inline static autoPtr<eConst> New(const dictionary& dict);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "eConst<"
              + word(EquationOfState::typeName()) + '>';
        }

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

        //- Absolute internal energy [J/kg]
        inline scalar Ea
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

        inline void operator+=(const eConst&);


    // Friend operators

        friend eConst operator+ <EquationOfState>
        (
            const eConst&,
            const eConst&
        );

        friend eConst operator* <EquationOfState>
        (
            const scalar,
            const eConst&
        );

        friend eConst operator== <EquationOfState>
        (
            const eConst&,
            const eConst&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include"eConstI.H"

#ifdef NoRepository
    #include "eConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
