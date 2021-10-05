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
    Foam::BirchMurnaghan2

Description
    Second order BirchMurnaghan equation of state for a solid

    References:
    \verbatim
        Zheng, Z., Zhao, J. (2016).
        Unreacted equation of states of typical energetic material under static
        compression: A review
        Chin. Phys. B, 25 (7), 1-11.
    \endverbatim

SourceFiles
    BirchMurnaghan2.C

\*---------------------------------------------------------------------------*/

#ifndef BirchMurnaghan2Blast_H
#define BirchMurnaghan2Blast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class BirchMurnaghan2;

template<class Specie>
inline BirchMurnaghan2<Specie> operator+
(
    const BirchMurnaghan2<Specie>&,
    const BirchMurnaghan2<Specie>&
);

template<class Specie>
inline BirchMurnaghan2<Specie> operator*
(
    const scalar,
    const BirchMurnaghan2<Specie>&
);

template<class Specie>
inline BirchMurnaghan2<Specie> operator==
(
    const BirchMurnaghan2<Specie>&,
    const BirchMurnaghan2<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const BirchMurnaghan2<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class BirchMurnaghan Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class BirchMurnaghan2
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- Reference pressure
    scalar pRef_;

    //- Bulk modulus
    scalar K0_;

    //- Derivative of bulk modulus w.r.t. pressure
    scalar K0Prime_;

    //- Mie Grunseisen coefficient
    scalar Gamma_;


public:

    // Constructors

        //- Construct from components
        inline BirchMurnaghan2
        (
            const Specie& sp,
            const scalar rho0,
            const scalar pRef,
            const scalar K0,
            const scalar K0Prime,
            const scalar Gamma
        );

        //- Construct from dictionary
        BirchMurnaghan2(const dictionary& dict);

        //- Construct as named copy
        inline BirchMurnaghan2
        (
            const word& name,
            const BirchMurnaghan2&
        );

        //- Construct and return a clone
        inline autoPtr<BirchMurnaghan2> clone() const;

        // Selector from dictionary
        inline static autoPtr<BirchMurnaghan2> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "BirchMurnaghan2<"
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

        //- Gamma function
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Speed of sound
        inline scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return derivative of pressure w.r.t. internal energy
        inline scalar dpde
        (
            const scalar rho,
            const scalar e,
            const scalar T
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

        inline void operator+=(const BirchMurnaghan2&);
        inline void operator*=(const scalar);


    // Friend operators

        friend BirchMurnaghan2 operator+ <Specie>
        (
            const BirchMurnaghan2&,
            const BirchMurnaghan2&
        );

        friend BirchMurnaghan2 operator* <Specie>
        (
            const scalar s,
            const BirchMurnaghan2&
        );

        friend BirchMurnaghan2 operator== <Specie>
        (
            const BirchMurnaghan2&,
            const BirchMurnaghan2&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const BirchMurnaghan2&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BirchMurnaghan2I.H"

#ifdef NoRepository
    #include "BirchMurnaghan2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
