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
    Foam::AbelNobel

Description
    Abel-Nobel equation of state for ballistics

    References:
    \verbatim
        Johnston, I.A., (2005).
        The Noble-Abel Equation of State: Thermodynamic Derivations for
        Ballistics Modeling (No. DSTO-TN-0670).
        Defense Science and Technology Organisation.
    \endverbatim

SourceFiles
    AbelNobel.C

\*---------------------------------------------------------------------------*/

#ifndef AbelNobelBlast_H
#define AbelNobelBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class AbelNobel;

template<class Specie>
inline AbelNobel<Specie> operator+
(
    const AbelNobel<Specie>&,
    const AbelNobel<Specie>&
);

template<class Specie>
inline AbelNobel<Specie> operator*
(
    const scalar,
    const AbelNobel<Specie>&
);

template<class Specie>
inline AbelNobel<Specie> operator==
(
    const AbelNobel<Specie>&,
    const AbelNobel<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const AbelNobel<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class AbelNobel Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class AbelNobel
:
    public Specie
{
// Private data

    //- Coefficients
    scalar b_;


public:

    // Constructors

        //- Construct from components
        inline AbelNobel
        (
            const Specie& sp,
            const scalar b
        );

        //- Construct from dictionary
        AbelNobel(const dictionary& dict);

        //- Construct as named copy
        inline AbelNobel
        (
            const word& name,
            const AbelNobel&
        );

        //- Construct and return a clone
        inline autoPtr<AbelNobel> clone() const;

        // Selector from dictionary
        inline static autoPtr<AbelNobel> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "AbelNobel<"
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

        //- Return the MieGruniesen coefficient
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return the speed of sound
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
            const scalar T
        ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const AbelNobel&);
        inline void operator*=(const scalar);


    // Friend operators

        friend AbelNobel operator+ <Specie>
        (
            const AbelNobel&,
            const AbelNobel&
        );

        friend AbelNobel operator* <Specie>
        (
            const scalar s,
            const AbelNobel&
        );

        friend AbelNobel operator== <Specie>
        (
            const AbelNobel&,
            const AbelNobel&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const AbelNobel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AbelNobelI.H"

#ifdef NoRepository
    #include "AbelNobel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
