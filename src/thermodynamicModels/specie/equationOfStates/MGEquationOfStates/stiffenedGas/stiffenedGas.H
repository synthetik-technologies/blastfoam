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
    Foam::stiffenedGas

Description
    Equation of state for a weakly compressible material

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    stiffenedGas.C

\*---------------------------------------------------------------------------*/

#ifndef stiffenedGasBlast_H
#define stiffenedGasBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class stiffenedGas;

template<class Specie>
inline stiffenedGas<Specie> operator+
(
    const stiffenedGas<Specie>&,
    const stiffenedGas<Specie>&
);

template<class Specie>
inline stiffenedGas<Specie> operator*
(
    const scalar,
    const stiffenedGas<Specie>&
);

template<class Specie>
inline stiffenedGas<Specie> operator==
(
    const stiffenedGas<Specie>&,
    const stiffenedGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const stiffenedGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class stiffenedGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class stiffenedGas
:
    public Specie
{
// Private data

        //- Model coefficients
        scalar a_;

        //- Specific heat ratio
        scalar gamma_;

        //- Minimum density used for internal energy contribution
        scalar minRho_;


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
        inline stiffenedGas
        (
            const Specie& sp,
            const scalar a,
            const scalar gamma
        );

        //- Construct from dictionary
        stiffenedGas(const dictionary& dict);

        //- Construct as named copy
        inline stiffenedGas
        (
            const word& name,
            const stiffenedGas&
        );

        //- Construct and return a clone
        inline autoPtr<stiffenedGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<stiffenedGas> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "stiffenedGas<"
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

        //- Return density given temperature and pressure
        inline scalar rho(const scalar p, const scalar T) const;

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

        //- ENtropy of mixture [J/kg/K]
        inline scalar S
        (
            const scalar p,
            const scalar T
        ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const stiffenedGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend stiffenedGas operator+ <Specie>
        (
            const stiffenedGas&,
            const stiffenedGas&
        );

        friend stiffenedGas operator* <Specie>
        (
            const scalar s,
            const stiffenedGas&
        );

        friend stiffenedGas operator== <Specie>
        (
            const stiffenedGas&,
            const stiffenedGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const stiffenedGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stiffenedGasI.H"

#ifdef NoRepository
    #include "stiffenedGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
