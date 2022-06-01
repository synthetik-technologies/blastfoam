/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::Murnaghan

Description
    Murnaghan equation of state for a solid

    References:
    \verbatim
        Souers, P.C., Anderson, S., Mercer, J., McGuire, E., Vitello, P.
        (2000).
        JWL++: A Simple Reactive Flow Code Package for Detonation
        Propellants, Explosives, Pyrotechnics, 25, 54-58.
    \endverbatim

SourceFiles
    Murnaghan.C

\*---------------------------------------------------------------------------*/

#ifndef MurnaghanBlast_H
#define MurnaghanBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class Murnaghan;

template<class Specie>
inline Murnaghan<Specie> operator+
(
    const Murnaghan<Specie>&,
    const Murnaghan<Specie>&
);

template<class Specie>
inline Murnaghan<Specie> operator*
(
    const scalar,
    const Murnaghan<Specie>&
);

template<class Specie>
inline Murnaghan<Specie> operator==
(
    const Murnaghan<Specie>&,
    const Murnaghan<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const Murnaghan<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class Murnaghan Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Murnaghan
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- Reference pressure
    scalar pRef_;

    //- Coefficients
    scalar K0_;
    scalar kappa_;

    scalar K0Prime_;
    scalar n_;

    //- Mie Grunsien coefficient
    scalar Gamma_;


public:

    // Constructors

        //- Construct from components
        inline Murnaghan
        (
            const Specie& sp,
            const scalar rho0,
            const scalar pRef,
            const scalar K0,
            const scalar kappa,
            const scalar K0Prime,
            const scalar n,
            const scalar Gamma
        );

        //- Construct from dictionary
        Murnaghan(const dictionary& dict);

        //- Construct as named copy
        inline Murnaghan
        (
            const word& name,
            const Murnaghan&
        );

        //- Construct and return a clone
        inline autoPtr<Murnaghan> clone() const;

        // Selector from dictionary
        inline static autoPtr<Murnaghan> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "Murnaghan<"
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
            return true;
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

        //- Gamma function
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return non-ideal speed of sound
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

        inline void operator+=(const Murnaghan&);
        inline void operator*=(const scalar);


    // Friend operators

        friend Murnaghan operator+ <Specie>
        (
            const Murnaghan&,
            const Murnaghan&
        );

        friend Murnaghan operator* <Specie>
        (
            const scalar s,
            const Murnaghan&
        );

        friend Murnaghan operator== <Specie>
        (
            const Murnaghan&,
            const Murnaghan&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const Murnaghan&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MurnaghanI.H"

#ifdef NoRepository
    #include "Murnaghan.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
