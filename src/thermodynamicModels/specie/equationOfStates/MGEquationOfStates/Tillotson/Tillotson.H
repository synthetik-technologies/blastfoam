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
    Tillotson equation of state

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
#include "lookupTable1D.H"

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


    // Cavitation parameters

        //- Is cavitation used
        Switch cavitation_;

        //- Cavitation pressure
        scalar pCav_;

        //- Cavitation density
        scalar rhoCav_;

        //- Ideal gas constant for vapor
        scalar k_;

        //- Stored cavitation energies
        lookupTable1D EcTable_;


public:

    //- Construct from components
        inline Tillotson
        (
            const Specie& sp,
            const scalar p0,
            const scalar pCav,
            const scalar rho0,
            const scalar rhoCav,
            const scalar e0,
            const scalar omega,
            const scalar k,
            const scalar A,
            const scalar B,
            const scalar C
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

        //- Is this a solid equation of state
        static bool solid()
        {
            return false;
        }

        //- Is temperature used
        bool temperatureBased() const
        {
            return false;
        }

        //- Return pressure
        scalar p
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const bool limit = true
        ) const;

        //- Return the speed of sound
        scalar cSqr
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const scalar& Cv
        ) const;


        //- Return deviation from ideal gas pressure
        scalar Pi
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Return non-ideal speed of sound
        scalar delta
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e
        ) const;

        //- Gamma function
        scalar Gamma
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const scalar& Cv
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. internal energy
        scalar dpde
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar& p,
            const scalar& rho,
            const scalar& T
        ) const;


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
