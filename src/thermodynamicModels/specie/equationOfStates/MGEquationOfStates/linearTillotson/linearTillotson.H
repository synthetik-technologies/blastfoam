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
    Foam::linearTillotson

Description
    linearTillotson equation of state

    References:
    \verbatim
        Wardlaw, A., Jr., McKnown, R., Luton, A.  (2011).
        Coupled Hydrocode Prediction of Underwater Explosion Damage.
        Warhead Performance and Lethality Division, Code 420
    \endverbatim

SourceFiles
    linearTillotson.C

\*---------------------------------------------------------------------------*/

#ifndef linearTillotson_H
#define linearTillotson_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "scalarLookupTable1D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class linearTillotson;

template<class Specie>
inline linearTillotson<Specie> operator+
(
    const linearTillotson<Specie>&,
    const linearTillotson<Specie>&
);

template<class Specie>
inline linearTillotson<Specie> operator*
(
    const scalar,
    const linearTillotson<Specie>&
);

template<class Specie>
inline linearTillotson<Specie> operator==
(
    const linearTillotson<Specie>&,
    const linearTillotson<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const linearTillotson<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class linearTillotson Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class linearTillotson
:
    public Specie
{
// Private data

    //-  Simplified Tillotson parameters

        //- Reference pressure
        scalar p0_;

        //- Reference density
        scalar rho0_;

        //- Reference internal energy
        scalar e0_;

        //- Coefficients
        scalar omega_;
        scalar A_;
        scalar B_;
        scalar C_;

        //- Cavitation pressure
        scalar pCav_;

        //- Stored cavitation energies
        scalarLookupTable1D EcTable_;


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
        inline linearTillotson
        (
            const Specie& sp,
            const scalar p0,
            const scalar rho0,
            const scalar e0,
            const scalar omega,
            const scalar A,
            const scalar B,
            const scalar C,
            const scalar pCav
        );

        //- Construct from dictionary
        linearTillotson(const dictionary& dict);

        //- Construct as named copy
        inline linearTillotson
        (
            const word& name,
            const linearTillotson&
        );

        //- Construct and return a clone
        inline autoPtr<linearTillotson> clone() const;

        // Selector from dictionary
        inline static autoPtr<linearTillotson> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "linearTillotson<"
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

        inline void operator+=(const linearTillotson&);
        inline void operator*=(const scalar);


    // Friend operators

        friend linearTillotson operator+ <Specie>
        (
            const linearTillotson&,
            const linearTillotson&
        );

        friend linearTillotson operator* <Specie>
        (
            const scalar s,
            const linearTillotson&
        );

        friend linearTillotson operator== <Specie>
        (
            const linearTillotson&,
            const linearTillotson&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const linearTillotson&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "linearTillotsonI.H"

#ifdef NoRepository
    #include "linearTillotson.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
