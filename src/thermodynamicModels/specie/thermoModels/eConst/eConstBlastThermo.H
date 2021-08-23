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
    Foam::eConstThermo

Description
    Thermodynamics model based with constant properties

SourceFiles
    eConstBlastThermoI.H
    eConstBlastThermo.C
\*---------------------------------------------------------------------------*/

#ifndef eConstThermoBlast_H
#define eConstThermoBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModel.H"
#include "CpCvModel.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eConstThermo;

template<class EquationOfState>
inline eConstThermo<EquationOfState> operator+
(
    const eConstThermo<EquationOfState>&,
    const eConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline eConstThermo<EquationOfState> operator*
(
    const scalar,
    const eConstThermo<EquationOfState>&
);

template<class EquationOfState>
inline eConstThermo<EquationOfState> operator==
(
    const eConstThermo<EquationOfState>&,
    const eConstThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eConstThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class eConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eConstThermo
:
    public EquationOfState,
    public CpCvModel
{
    // Private data

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Flame temperature [K]
        scalar flameT_;


        //- Construct from components
        inline eConstThermo
        (
            const EquationOfState& st,
            const scalar cv,
            const scalar hf,
            const scalar ft
        );


public:

    // Constructors

        //- Construct from dictionary
        eConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline eConstThermo(const word&, const eConstThermo&);

        //- Construct as copy
        inline eConstThermo(const eConstThermo&);

        //- Construct and return a clone
        inline autoPtr<eConstThermo> clone() const;

        // Selector from dictionary
        inline static autoPtr<eConstThermo> New(const dictionary& dict);


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "eConst<"
              + word(EquationOfState::typeName()) + '>';
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

        inline void operator=(const eConstThermo&);
        inline void operator+=(const eConstThermo&);


    // Friend operators

        friend eConstThermo operator+ <EquationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );

        friend eConstThermo operator* <EquationOfState>
        (
            const scalar,
            const eConstThermo&
        );

        friend eConstThermo operator== <EquationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include"eConstBlastThermoI.H"

#ifdef NoRepository
    #include "eConstBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
