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
    Foam::solidJWL

Description
    Copy of JWL equation of state used for a blended JWL detonating material

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    solidJWL.C

\*---------------------------------------------------------------------------*/

#ifndef solidJWLBlast_H
#define solidJWLBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class solidJWL;

template<class Specie>
inline solidJWL<Specie> operator+
(
    const solidJWL<Specie>&,
    const solidJWL<Specie>&
);

template<class Specie>
inline solidJWL<Specie> operator*
(
    const scalar,
    const solidJWL<Specie>&
);

template<class Specie>
inline solidJWL<Specie> operator==
(
    const solidJWL<Specie>&,
    const solidJWL<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const solidJWL<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class solidJWL Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class solidJWL
:
    public Specie
{
// Private data

        //- Reference density
        scalar rho0_;

        //- Reference Gamma
        scalar omega_;

        //- Model coefficients
        scalar A_;
        scalar B_;
        scalar R1_;
        scalar R2_;

        //- Reference pressure
        scalar pRef_;

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
        inline solidJWL
        (
            const Specie& sp,
            const scalar rho0,
            const scalar omega,
            const scalar A,
            const scalar B,
            const scalar R1,
            const scalar R2,
            const scalar pRef,
            const scalar e0
        );

        //- Construct from dictionary
        solidJWL(const dictionary& dict);

        //- Construct as named copy
        inline solidJWL
        (
            const word& name,
            const solidJWL&
        );

        //- Construct and return a clone
        inline autoPtr<solidJWL> clone() const;

        // Selector from dictionary
        inline static autoPtr<solidJWL> New
        (
            const dictionary& dict
        );

    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "solidJWL<"
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

        inline void operator+=(const solidJWL&);
        inline void operator*=(const scalar);


    // Friend operators

        friend solidJWL operator+ <Specie>
        (
            const solidJWL&,
            const solidJWL&
        );

        friend solidJWL operator* <Specie>
        (
            const scalar s,
            const solidJWL&
        );

        friend solidJWL operator== <Specie>
        (
            const solidJWL&,
            const solidJWL&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const solidJWL&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidJWLI.H"

#ifdef NoRepository
    #include "solidJWL.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
