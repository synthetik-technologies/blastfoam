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
    Foam::JWL

Description
    Equation of state for an energetic material

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    JWL.C

\*---------------------------------------------------------------------------*/

#ifndef JWLBlast_H
#define JWLBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class JWL;

template<class Specie>
inline JWL<Specie> operator+
(
    const JWL<Specie>&,
    const JWL<Specie>&
);

template<class Specie>
inline JWL<Specie> operator*
(
    const scalar,
    const JWL<Specie>&
);

template<class Specie>
inline JWL<Specie> operator==
(
    const JWL<Specie>&,
    const JWL<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const JWL<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class JWL Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class JWL
:
    public Specie
{
// Private data

        //- Reference density
        scalar rho0_;

        //- Cut off density (exponents are no longer considered)
        scalar rhoCutOff_;

        //- Reference Gamma
        scalar omega_;

        //- Ideal gas specific heat ratio
        scalar gammaIdeal_;

        //- Model coefficients
        scalar A_;
        scalar B_;
        scalar R1_;
        scalar R2_;

        //- Internal energy that gives 0 pressure at reference density
        scalar e0_;


    // Private member functions

        //- Return deviation from ideal gas pressure
        inline scalar Pi
        (
            const scalar rho,
            const scalar e
        ) const;

        //- Return non-ideal speed of sound
        inline scalar delta
        (
            const scalar p,
            const scalar rho,
            const scalar e
        ) const;


public:

    // Constructors

        //- Construct from components
        inline JWL
        (
            const Specie& sp,
            const scalar rho0,
            const scalar rhoCutOff,
            const scalar omega,
            const scalar gammaIdeal,
            const scalar A,
            const scalar B,
            const scalar R1,
            const scalar R2,
            const scalar e0
        );

        //- Construct from dictionary
        JWL(const dictionary& dict);

        //- Construct as named copy
        inline JWL
        (
            const word& name,
            const JWL&
        );

        //- Construct and return a clone
        inline autoPtr<JWL> clone() const;

        // Selector from dictionary
        inline static autoPtr<JWL> New
        (
            const dictionary& dict
        );

    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "JWL<"
            + word(Specie::typeName_()) + '>';
        }

        //- Is temperature used
        static inline bool temperatureBased()
        {
            return false;
        }

        //- Is this a solid equation of state
        static inline bool solid()
        {
            return false;
        }

        //- Return pressure
        inline scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const bool limit = true
        ) const;

        //- Return density given temperature and pressure
        inline scalar rho(const scalar p, const scalar T) const;

        //- Return the MieGruniesen coefficient
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar Cv
        ) const;

        //- Return the speed of sound
        inline scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar Cv
        ) const;

        //- Return derivative of pressure w.r.t density
        inline scalar dpdv
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

        //- Return derivative of pressure w.r.t. Temperature
        inline scalar dpdT
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

        inline void operator+=(const JWL&);
        inline void operator*=(const scalar);


    // Friend operators

        friend JWL operator+ <Specie>
        (
            const JWL&,
            const JWL&
        );

        friend JWL operator* <Specie>
        (
            const scalar s,
            const JWL&
        );

        friend JWL operator== <Specie>
        (
            const JWL&,
            const JWL&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const JWL&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "JWLI.H"

#ifdef NoRepository
    #include "JWL.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
