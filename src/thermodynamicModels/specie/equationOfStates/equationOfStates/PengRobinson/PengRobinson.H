
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
        Vc           | Critical volume [m^3/kmol]
        Pc           | Critical pressure [Pa]
        omega        | Acentric factor [-]
    \endtable

    Example specification of the PengRobinsonGas equation of state for methane:

    \verbatim
    equationOfState
    {
        Tc        190.55;
        Vc        0.0285;
        Pc        4.595e6;
        omega     0.0115;
    }
    \endverbatim

SourceFiles
    PengRobinsonI.H
    PengRobinson.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonBlast_H
#define PengRobinsonBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class PengRobinson;

template<class Specie>
inline PengRobinson<Specie> operator+
(
    const PengRobinson<Specie>&,
    const PengRobinson<Specie>&
);

template<class Specie>
inline PengRobinson<Specie> operator*
(
    const scalar,
    const PengRobinson<Specie>&
);

template<class Specie>
inline PengRobinson<Specie> operator==
(
    const PengRobinson<Specie>&,
    const PengRobinson<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const PengRobinson<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class PengRobinson Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class PengRobinson
:
    public Specie
{
    // Private Data

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Saved parameters
        scalar a_;
        scalar b_;
        scalar kappa_;


public:

    // Constructors

        //- Construct from components
        inline PengRobinson
        (
            const Specie& sp,
            const scalar Tc,
            const scalar Vc,
            const scalar Pc,
            const scalar omega,
            const scalar Zc
        );

        //- Construct from dictionary
        PengRobinson(const dictionary& dict);

        //- Construct as named copy
        inline PengRobinson
        (
            const word& name,
            const PengRobinson&
        );

        //- Construct and return a clone
        inline autoPtr<PengRobinson> clone() const;

        // Selector from dictionary
        inline static autoPtr<PengRobinson> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "PengRobinson<"
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

        //- Return compressibility
        inline scalar Z(const scalar p, const scalar T) const;

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

        inline void operator+=(const PengRobinson&);
        inline void operator*=(const scalar);


    // Friend operators

        friend PengRobinson operator+ <Specie>
        (
            const PengRobinson&,
            const PengRobinson&
        );

        friend PengRobinson operator* <Specie>
        (
            const scalar s,
            const PengRobinson&
        );

        friend PengRobinson operator== <Specie>
        (
            const PengRobinson&,
            const PengRobinson&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const PengRobinson&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonI.H"

#ifdef NoRepository
    #include "PengRobinson.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
