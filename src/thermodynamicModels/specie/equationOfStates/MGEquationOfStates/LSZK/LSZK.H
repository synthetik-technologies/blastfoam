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
    Foam::LSZK

Description
    Simple equation of state for an energetic material

    References:
    \verbatim
        Needham, C. E.  (2018).
        Shock Wave and High Pressure Phenomena.
        Springer
    \endverbatim

SourceFiles
    LSZK.C

\*---------------------------------------------------------------------------*/

#ifndef LSZK_H
#define LSZK_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class LSZK;

template<class Specie>
inline LSZK<Specie> operator+
(
    const LSZK<Specie>&,
    const LSZK<Specie>&
);

template<class Specie>
inline LSZK<Specie> operator*
(
    const scalar,
    const LSZK<Specie>&
);

template<class Specie>
inline LSZK<Specie> operator==
(
    const LSZK<Specie>&,
    const LSZK<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const LSZK<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class LSZK Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class LSZK
:
    public Specie
{
// Private data

        //- Model coefficients
        scalar a_;
        scalar b_;

        //- Specific heat ratio
        scalar gamma_;


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
        inline LSZK
        (
            const Specie& sp,
            const scalar a,
            const scalar b,
            const scalar gamma
        );

        //- Construct from dictionary
        LSZK(const dictionary& dict);

        //- Construct as named copy
        inline LSZK
        (
            const word& name,
            const LSZK&
        );

        //- Construct and return a clone
        inline autoPtr<LSZK> clone() const;

        // Selector from dictionary
        inline static autoPtr<LSZK> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "LSZK<"
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

        inline void operator+=(const LSZK&);
        inline void operator*=(const scalar);


    // Friend operators

        friend LSZK operator+ <Specie>
        (
            const LSZK&,
            const LSZK&
        );

        friend LSZK operator* <Specie>
        (
            const scalar s,
            const LSZK&
        );

        friend LSZK operator== <Specie>
        (
            const LSZK&,
            const LSZK&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const LSZK&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LSZKI.H"

#ifdef NoRepository
    #include "LSZK.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
