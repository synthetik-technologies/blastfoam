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
#include "MGEquationOfState.H"
#include "lookupTable1D.H"

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
    public Specie,
    public MGEquationOfState
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
        lookupTable1D EcTable_;


public:

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
        virtual scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const bool limit = true
        ) const;

        //- Return deviation from ideal gas pressure
        virtual scalar Pi
        (
            const scalar rho,
            const scalar e
        ) const;

        //- Return non-ideal speed of sound
        virtual scalar delta
        (
            const scalar p,
            const scalar rho,
            const scalar e
        ) const;

        //- Gamma function
        virtual scalar GammaMG
        (
            const scalar rho,
            const scalar e
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t. internal energy
        scalar dpde
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
        virtual scalar CpMCv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const
        {
            return MGEquationOfState::CpMCv(this->cpcv(), rho, e, T);
        }

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar p,
            const scalar rho,
            const scalar T
        ) const;


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
