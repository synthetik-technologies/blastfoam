
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2022-07-19 Jeff Heylmun:    Modified class for a density based thermodynamic
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
    Foam::PengRobinsonGas

Description
    PengRobinsonGas cubic equation of state
    Coefficient mixing is very inaccurate and not supported, so this equation of
    state can only be applied to mixtures by using a thermophysical properties
    model which mixes specie property value.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Pc           | Critical pressure [Pa]
        a            | Atteaction factor [-]
        b            | Covolume [m^3]
    \endtable

    Example specification of the PengRobinsonGas equation of state for methane:

    \verbatim
    equationOfState
    {
        Tc        190.55;
        Pc        4.595e6;
    }
    \endverbatim

SourceFiles
    VanderWaalsI.H
    VanderWaals.C

\*---------------------------------------------------------------------------*/

#ifndef VanderWaals_H
#define VanderWaals_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class VanderWaals;

template<class Specie>
inline VanderWaals<Specie> operator+
(
    const VanderWaals<Specie>&,
    const VanderWaals<Specie>&
);

template<class Specie>
inline VanderWaals<Specie> operator*
(
    const scalar,
    const VanderWaals<Specie>&
);

template<class Specie>
inline VanderWaals<Specie> operator==
(
    const VanderWaals<Specie>&,
    const VanderWaals<Specie>&
);

template<class Specie>
inline Ostream& operator<<
(
    Ostream&,
    const VanderWaals<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class VanderWaals Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class VanderWaals
:
    public Specie
{
    // Private Data

        //- Attraction [K]
        scalar a_;

        //- Covolume [m^3]
        scalar b_;

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

public:

    // Constructors

        //- Construct from components
        inline VanderWaals
        (
            const Specie& sp,
            const scalar Tc,
            const scalar Pc
        );

        //- Construct from dictionary
        VanderWaals(const dictionary& dict);

        //- Construct as named copy
        inline VanderWaals
        (
            const word& name,
            const VanderWaals&
        );

        //- Construct and return a clone
        inline autoPtr<VanderWaals> clone() const;

        // Selector from dictionary
        inline static autoPtr<VanderWaals> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "VanderWaals<"
               + word(Specie::typeName_()) + '>';
        }

        //- Is temperature used
        static bool temperatureBased()
        {
            return true;
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return false;
        }

        //- Return pressure optionally limited pressure
        inline scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const bool limit = true
        ) const;

        //- Return unlimited pressure
        inline scalar pRhoT
        (
            const scalar rho,
            const scalar T
        ) const;

        //- Return density given temperature and pressure
        inline scalar rho(const scalar p, const scalar T) const;

        //- Return the MieGruniesen coefficient
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return the speed of sound
        inline scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return derivative of pressure w.r.t density
        inline scalar dpdv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        inline scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t. internal energy
        inline scalar dpde
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Internal energy correction [J/kg]
        inline scalar E
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Enthalpy correction [J/kg]
        inline scalar H
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        inline scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        inline scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity difference [J/kg/K]
        inline scalar CpMCv
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar CpCv
        ) const;

        //- Entropy of mixture [J/kg/K]
        inline scalar S
        (
            const scalar p,
            const scalar T
        ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const VanderWaals&);
        inline void operator*=(const scalar);


    // Friend operators

        friend VanderWaals operator+ <Specie>
        (
            const VanderWaals&,
            const VanderWaals&
        );

        friend VanderWaals operator* <Specie>
        (
            const scalar s,
            const VanderWaals&
        );

        friend VanderWaals operator== <Specie>
        (
            const VanderWaals&,
            const VanderWaals&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const VanderWaals&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VanderWaalsI.H"

#ifdef NoRepository
    #include "VanderWaals.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
