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
    Foam::BWR

Description
    Bennedict-Webb-Rubin equation of state

    References:
    \verbatim
        Bennedict, M., Webb, G. B., Rubin, L.C. (1940).
        An Empirical Equation for Thermodynamic Properties of Light
        Hydrocarbons and Their Mixtures: I. Methane, Ethane, Propane, and
        n-Butane
        Journal of Chemical Physics, 8 (4), 334-345.
    \endverbatim

SourceFiles
    BWR.C

\*---------------------------------------------------------------------------*/

#ifndef BWRBlast_H
#define BWRBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class BWR;

template<class Specie>
inline BWR<Specie> operator+
(
    const BWR<Specie>&,
    const BWR<Specie>&
);

template<class Specie>
inline BWR<Specie> operator*
(
    const scalar,
    const BWR<Specie>&
);

template<class Specie>
inline BWR<Specie> operator==
(
    const BWR<Specie>&,
    const BWR<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const BWR<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class BWR Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class BWR
:
    public Specie
{
// Private data

    //- Coefficients
    scalar A0_;
    scalar B0_;
    scalar C0_;
    scalar a_;
    scalar b_;
    scalar c_;
    scalar alpha_;
    scalar rhoc_;
    scalar gamma_;


public:

    // Constructors

        //- Construct from components
        inline BWR
        (
            const Specie& sp,
            const scalar A0,
            const scalar B0,
            const scalar C0,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar alpha,
            const scalar rhoc,
            const scalar gamma
        );

        //- Construct from dictionary
        BWR(const dictionary& dict);

        //- Construct as named copy
        inline BWR
        (
            const word& name,
            const BWR&
        );

        //- Construct and return a clone
        inline autoPtr<BWR> clone() const;

        // Selector from dictionary
        inline static autoPtr<BWR> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "BWR<"
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

        //- Heat capacity at constant volume of mixture [J/kg/K]
        inline scalar S
        (
            const scalar p,
            const scalar rho,
            const scalar T
        ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const BWR&);
        inline void operator*=(const scalar);


    // Friend operators

        friend BWR operator+ <Specie>
        (
            const BWR&,
            const BWR&
        );

        friend BWR operator* <Specie>
        (
            const scalar s,
            const BWR&
        );

        friend BWR operator== <Specie>
        (
            const BWR&,
            const BWR&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const BWR&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BWRI.H"

#ifdef NoRepository
    #include "BWR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
