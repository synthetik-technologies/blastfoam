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
    Foam::JWLC

Description
    C form of the JWL equation of state (no internal energy dependence)

    References:
    \verbatim
        Souers, P.C., Anderson, S., Mercer, J., McGuire, E., Vitello, P.
        (2000).
        JWL++: A Simple Reactive Flow Code Package for Detonation
        Propellants, Explosives, Pyrotechnics, 25, 54-58.
    \endverbatim

SourceFiles
    JWLC.C

\*---------------------------------------------------------------------------*/

#ifndef JWLCBlast_H
#define JWLCBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "equationOfState.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class JWLC;

template<class Specie>
inline JWLC<Specie> operator+
(
    const JWLC<Specie>&,
    const JWLC<Specie>&
);

template<class Specie>
inline JWLC<Specie> operator*
(
    const scalar,
    const JWLC<Specie>&
);

template<class Specie>
inline JWLC<Specie> operator==
(
    const JWLC<Specie>&,
    const JWLC<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const JWLC<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class JWLC Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class JWLC
:
    public Specie,
    public equationOfState
{
// Private data

        //- Reference density
        scalar rho0_;

        //- Mie Grunsien coefficient
        scalar omega_;

        //- Model coefficients
        scalar A_;
        scalar B_;
        scalar C_;
        scalar R1_;
        scalar R2_;



public:

    // Constructors

        //- Construct from components
        inline JWLC
        (
            const Specie& sp,
            const scalar rho0,
            const scalar omega,
            const scalar A,
            const scalar B,
            const scalar C,
            const scalar R1,
            const scalar R2
        );

        //- Construct from dictionary
        JWLC(const dictionary& dict);

        //- Construct as named copy
        inline JWLC
        (
            const word& name,
            const JWLC&
        );

        //- Construct and return a clone
        inline autoPtr<JWLC> clone() const;

        // Selector from dictionary
        inline static autoPtr<JWLC> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "JWLC<"
            + word(Specie::typeName_()) + '>';
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return false;
        }

        //- Return deviation from ideal gas pressure
        virtual scalar pRhoT
        (
            const scalar rho,
            const scalar T
        ) const;

        //- Return non-ideal speed of sound
        virtual scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Gamma function
        virtual scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar p,
            const scalar rho,
            const scalar T
        ) const;


    // Member operators

        inline void operator+=(const JWLC&);
        inline void operator*=(const scalar);


    // Friend operators

        friend JWLC operator+ <Specie>
        (
            const JWLC&,
            const JWLC&
        );

        friend JWLC operator* <Specie>
        (
            const scalar s,
            const JWLC&
        );

        friend JWLC operator== <Specie>
        (
            const JWLC&,
            const JWLC&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "JWLCI.H"

#ifdef NoRepository
    #include "JWLC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
