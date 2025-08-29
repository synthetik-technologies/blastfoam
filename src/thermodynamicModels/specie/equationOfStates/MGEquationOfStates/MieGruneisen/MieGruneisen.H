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
    Foam::MieGruneisen

Description
    Equation of state for an ideal gas using the Mie Grunsien form

SourceFiles
    MieGruneisen.C

\*---------------------------------------------------------------------------*/

#ifndef MieGruneisenBlast_H
#define MieGruneisenBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class MieGruneisen;

template<class Specie>
inline MieGruneisen<Specie> operator+
(
    const MieGruneisen<Specie>&,
    const MieGruneisen<Specie>&
);

template<class Specie>
inline MieGruneisen<Specie> operator*
(
    const scalar,
    const MieGruneisen<Specie>&
);

template<class Specie>
inline MieGruneisen<Specie> operator==
(
    const MieGruneisen<Specie>&,
    const MieGruneisen<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const MieGruneisen<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class MieGruneisen Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class MieGruneisen
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- Pressure offset
    scalar pInf_;

    //- Reference speed of sound
    scalar c0_;

    //- Grunsien parameter
    scalar Gamma_;

    // Private member functions

        //- Return deviation from ideal gas pressure
        inline scalar Pi(const scalar rho) const;

public:

    // Constructors

        //- Construct from components
        inline MieGruneisen
        (
            const Specie& sp,
            const scalar rho0,
            const scalar pInf,
            const scalar c0,
            const scalar Gamma
        );

        //- Construct from dictionary
        MieGruneisen(const dictionary& dict);

        //- Construct as named copy
        inline MieGruneisen
        (
            const word& name,
            const MieGruneisen&
        );

        //- Construct and return a clone
        inline autoPtr<MieGruneisen> clone() const;

        // Selector from dictionary
        inline static autoPtr<MieGruneisen> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "MieGruneisen<"
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

        //- Return energy given pressure and density
        inline scalar ePRho(const scalar p, const scalar rho) const;

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

        //- Return entropy contribution to the integral of Cp/T [J/kg/K]
        inline scalar Sp(const scalar p, const scalar T) const;

        //- Return entropy contribution to the integral of Cv/T [J/kg/K]
        inline scalar Sv(const scalar p, const scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const MieGruneisen&);
        inline void operator*=(const scalar);


    // Friend operators

        friend MieGruneisen operator+ <Specie>
        (
            const MieGruneisen&,
            const MieGruneisen&
        );

        friend MieGruneisen operator* <Specie>
        (
            const scalar s,
            const MieGruneisen&
        );

        friend MieGruneisen operator== <Specie>
        (
            const MieGruneisen&,
            const MieGruneisen&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const MieGruneisen&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MieGruneisenI.H"

#ifdef NoRepository
    #include "MieGruneisen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
