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
    Foam::tabulatedThermoEOS

Description
    Combination thermodynamic model and equation of state that uses
    tables to return state variables

SourceFiles
    tabulatedThermoEOSI.H
    tabulatedThermoEOS.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedThermoEOS_H
#define tabulatedThermoEOS_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTable2D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class tabulatedThermoEOS;

template<class Specie>
inline tabulatedThermoEOS<Specie> operator+
(
    const tabulatedThermoEOS<Specie>&,
    const tabulatedThermoEOS<Specie>&
);

template<class Specie>
inline tabulatedThermoEOS<Specie> operator*
(
    const scalar,
    const tabulatedThermoEOS<Specie>&
);

template<class Specie>
inline tabulatedThermoEOS<Specie> operator==
(
    const tabulatedThermoEOS<Specie>&,
    const tabulatedThermoEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const tabulatedThermoEOS<Specie>&
);

/*---------------------------------------------------------------------------*\
                    Class tabulatedThermoEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class tabulatedThermoEOS
:
    public Specie
{
// Private data

    //- Pressure lookup table
    lookupTable2D pTable_;

    //- Temperature lookup table
    lookupTable2D TTable_;


public:

    //- Construct from components
        inline tabulatedThermoEOS
        (
            const Specie& sp
        );

        //- Construct from dictionary
        tabulatedThermoEOS(const dictionary& dict);

        //- Construct as named copy
        inline tabulatedThermoEOS
        (
            const word& name,
            const tabulatedThermoEOS&
        );

        //- Construct and return a clone
        inline autoPtr<tabulatedThermoEOS> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedThermoEOS> New
        (
            const dictionary& dict
        );


    //- Return the instantiated type name
        static word typeName()
        {
            return
                "tabulated<tabulated<"
              + word(Specie::typeName_()) + ">>";
        }


    // Member Functions

        //- Is temperature used
        bool temperatureBased() const
        {
            return false;
        }

        // Equation of state functions

            //- Return pressure
            scalar p(const scalar& rho, const scalar& e, const scalar& T) const;

            //- Return derivative of pressure w.r.t. specific volume
            scalar dpdv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return derivative of pressure w.r.t. Temperature
            scalar dpdT
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return derivative of pressure w.r.t. internal energy
            scalar dpde
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

        // Thermodynamic functions

            //- Return Mie Gruniesen coefficient
            scalar Gamma
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return speed of sound
            scalar cSqr
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            scalar Cv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            scalar Cp
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Heat capacity difference [J/kg/K]
            scalar CpMCv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return specific heat ratio
            scalar CpByCv
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Sensible internal energy [J/kg]
            scalar Es
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Absolute internal energy [J/kg]
            scalar Ea
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Sensible enthalpy [J/kg]
            scalar Hs
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Absolute enthalpy [J/kg]
            scalar Ha
            (
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Heat of formation [J/kg]
            scalar Hf() const;

            //- Entropy [J/kg/K]
            scalar S
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Return temperature
            scalar TRhoE
            (
                const scalar& T0,
                const scalar& rho,
                const scalar& e
            ) const;

            //- Initialize internal energy
            scalar initializeEnergy
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

            //- Initialize internal energy
            scalar initializeRho
            (
                const scalar& p,
                const scalar& rho,
                const scalar& e,
                const scalar& T
            ) const;

        // Member operators

        inline void operator+=(const tabulatedThermoEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend tabulatedThermoEOS operator+ <Specie>
        (
            const tabulatedThermoEOS&,
            const tabulatedThermoEOS&
        );

        friend tabulatedThermoEOS operator* <Specie>
        (
            const scalar s,
            const tabulatedThermoEOS&
        );

        friend tabulatedThermoEOS operator== <Specie>
        (
            const tabulatedThermoEOS&,
            const tabulatedThermoEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedThermoEOSI.H"

#ifdef NoRepository
    #include "tabulatedThermoEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
