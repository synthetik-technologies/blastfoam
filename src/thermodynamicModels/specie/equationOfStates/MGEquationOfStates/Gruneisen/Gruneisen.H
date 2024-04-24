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
    Foam::Gruneisen

Description
    Equation of state for an ideal gas using the Mie Grunsien form

SourceFiles
    Gruneisen.C

\*---------------------------------------------------------------------------*/

#ifndef GruneisenBlast_H
#define GruneisenBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class Gruneisen;

template<class Specie>
inline Gruneisen<Specie> operator+
(
    const Gruneisen<Specie>&,
    const Gruneisen<Specie>&
);

template<class Specie>
inline Gruneisen<Specie> operator*
(
    const scalar,
    const Gruneisen<Specie>&
);

template<class Specie>
inline Gruneisen<Specie> operator==
(
    const Gruneisen<Specie>&,
    const Gruneisen<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const Gruneisen<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class Gruneisen Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Gruneisen
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- vp-vs intersept
    scalar C_;

    //- Coefficients
    scalar S1_;
    scalar S2_;
    scalar S3_;

    //- Gruniesen gamma
    scalar Gamma_;

    //- Gamma correction
    scalar a_;

    // Private member functions

        //- Return deviation from ideal gas pressure
        inline scalar Pi
        (
            const scalar rho,
            const scalar e
        ) const;

        //- Return derivative of deviation from ideal gas pressure
        inline scalar dPidRho
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

        //- Derivative of Gamma function
        inline scalar dGammaMGdRho
        (
            const scalar rho,
            const scalar e
        ) const;


public:

    // Constructors

        //- Construct from components
        inline Gruneisen
        (
            const Specie& sp,
            const scalar rho0,
            const scalar C,
            const scalar S1,
            const scalar S2,
            const scalar S3,
            const scalar Gamma,
            const scalar a
        );

        //- Construct from dictionary
        Gruneisen(const dictionary& dict);

        //- Construct as named copy
        inline Gruneisen
        (
            const word& name,
            const Gruneisen&
        );

        //- Construct and return a clone
        inline autoPtr<Gruneisen> clone() const;

        // Selector from dictionary
        inline static autoPtr<Gruneisen> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "Gruneisen<"
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

        //- Heat capacity at constant volume of mixture [J/kg/K]
        inline scalar S
        (
            const scalar p,
            const scalar T
        ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const Gruneisen&);
        inline void operator*=(const scalar);


    // Friend operators

        friend Gruneisen operator+ <Specie>
        (
            const Gruneisen&,
            const Gruneisen&
        );

        friend Gruneisen operator* <Specie>
        (
            const scalar s,
            const Gruneisen&
        );

        friend Gruneisen operator== <Specie>
        (
            const Gruneisen&,
            const Gruneisen&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const Gruneisen&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GruneisenI.H"

#ifdef NoRepository
    #include "Gruneisen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
