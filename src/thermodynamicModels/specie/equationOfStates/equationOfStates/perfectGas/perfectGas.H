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
    Foam::perfectGas

Description
    Ideal gas equation of state based on temperature

SourceFiles
    perfectGas.C
    perfectGasI.H

\*---------------------------------------------------------------------------*/

#ifndef perfectGasBlast_H
#define perfectGasBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "equationOfState.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class perfectGas;

template<class Specie>
inline perfectGas<Specie> operator+
(
    const perfectGas<Specie>&,
    const perfectGas<Specie>&
);

template<class Specie>
inline perfectGas<Specie> operator*
(
    const scalar,
    const perfectGas<Specie>&
);

template<class Specie>
inline perfectGas<Specie> operator==
(
    const perfectGas<Specie>&,
    const perfectGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const perfectGas<Specie>&
);

/*---------------------------------------------------------------------------*\
                           Class perfectGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class perfectGas
:
    public Specie,
    public equationOfState
{
public:

    // Constructors

        //- Construct from components
        inline perfectGas
        (
            const Specie& sp
        );

        //- Construct from dictionary
        perfectGas(const dictionary& dict);

        //- Construct as named copy
        inline perfectGas
        (
            const word& name,
            const perfectGas&
        );

        //- Construct and return a clone
        inline autoPtr<perfectGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<perfectGas> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "perfectGas<"
            + word(Specie::typeName_()) + '>';
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return false;
        }

        //- Return deviation from ideal gas pressure
        virtual scalar pRhoT
        (
            const scalar rho,
            const scalar T
        ) const;

        //- Return non-ideal speed of sound
        virtual scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Gamma function
        virtual scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar p,
            const scalar rho,
            const scalar T
        ) const;


    // Member operators

        inline void operator+=(const perfectGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend perfectGas operator+ <Specie>
        (
            const perfectGas&,
            const perfectGas&
        );

        friend perfectGas operator* <Specie>
        (
            const scalar s,
            const perfectGas&
        );

        friend perfectGas operator== <Specie>
        (
            const perfectGas&,
            const perfectGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "perfectGasI.H"

#ifdef NoRepository
    #include "perfectGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
