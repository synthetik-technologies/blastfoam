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
    Foam::tabulatedEOS

Description
    Tabulated equation of state

SourceFiles
    tabulatedEOS.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedEOS_H
#define tabulatedEOS_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "lookupTables2D.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class tabulatedEOS;

template<class Specie>
inline tabulatedEOS<Specie> operator+
(
    const tabulatedEOS<Specie>&,
    const tabulatedEOS<Specie>&
);

template<class Specie>
inline tabulatedEOS<Specie> operator*
(
    const scalar,
    const tabulatedEOS<Specie>&
);

template<class Specie>
inline tabulatedEOS<Specie> operator==
(
    const tabulatedEOS<Specie>&,
    const tabulatedEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const tabulatedEOS<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class tabulatedEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class tabulatedEOS
:
    public Specie
{
// Private data

    //- Pressure lookup table
    scalarLookupTable2D pTable_;


public:

    // Constructors

        //- Construct from dictionary
        tabulatedEOS(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedEOS> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedEOS> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "tabulatedEOS<"
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
            return false;
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

        //- Return density given temperature and pressure
        inline scalar rho(const scalar p, const scalar T) const;

        //- Gamma function
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return non-ideal speed of sound
        inline scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
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

        inline void operator+=(const tabulatedEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend tabulatedEOS operator+ <Specie>
        (
            const tabulatedEOS&,
            const tabulatedEOS&
        );

        friend tabulatedEOS operator* <Specie>
        (
            const scalar s,
            const tabulatedEOS&
        );

        friend tabulatedEOS operator== <Specie>
        (
            const tabulatedEOS&,
            const tabulatedEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const tabulatedEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedEOSI.H"

#ifdef NoRepository
    #include "tabulatedEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
