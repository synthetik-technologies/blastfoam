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
    Foam::Tillotson

Description
    Four state Tillotson equation of state for metallic materials

    References:
    \verbatim
        Wardlaw, A., Jr., McKnown, R., Luton, A.  (2011).
        Coupled Hydrocode Prediction of Underwater Explosion Damage.
        Warhead Performance and Lethality Division, Code 420
    \endverbatim

SourceFiles
    Tillotson.C

\*---------------------------------------------------------------------------*/

#ifndef Tillotson_H
#define Tillotson_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTables1D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class Tillotson;

template<class Specie>
inline Tillotson<Specie> operator+
(
    const Tillotson<Specie>&,
    const Tillotson<Specie>&
);

template<class Specie>
inline Tillotson<Specie> operator*
(
    const scalar,
    const Tillotson<Specie>&
);

template<class Specie>
inline Tillotson<Specie> operator==
(
    const Tillotson<Specie>&,
    const Tillotson<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const Tillotson<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class Tillotson Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Tillotson
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- Reference internal energy
    scalar e0_;

    //- Tillotson parameters
    scalar a_;
    scalar b_;
    scalar A_;
    scalar B_;

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

        //- Gamma function
        inline scalar GammaMG
        (
            const scalar rho,
            const scalar e
        ) const;


public:

    // Constructors

        //- Construct from components
        inline Tillotson
        (
            const Specie& sp,
            const scalar rho0,
            const scalar e0,
            const scalar a,
            const scalar b,
            const scalar A,
            const scalar B,
            const scalar pCav
        );

        //- Construct from dictionary
        Tillotson(const dictionary& dict);

        //- Construct as named copy
        inline Tillotson
        (
            const word& name,
            const Tillotson&
        );

        //- Construct and return a clone
        inline autoPtr<Tillotson> clone() const;

        // Selector from dictionary
        inline static autoPtr<Tillotson> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "Tillotson<"
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

        inline void operator+=(const Tillotson&);
        inline void operator*=(const scalar);


    // Friend operators

        friend Tillotson operator+ <Specie>
        (
            const Tillotson&,
            const Tillotson&
        );

        friend Tillotson operator* <Specie>
        (
            const scalar s,
            const Tillotson&
        );

        friend Tillotson operator== <Specie>
        (
            const Tillotson&,
            const Tillotson&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const Tillotson&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TillotsonI.H"

#ifdef NoRepository
    #include "Tillotson.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
