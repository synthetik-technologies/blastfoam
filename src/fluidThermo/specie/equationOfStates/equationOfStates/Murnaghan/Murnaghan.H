/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::Murnaghan

Description
    Murnaghan equation of state for a solid

    References:
    \verbatim
        Souers, P.C., Anderson, S., Mercer, J., McGuire, E., Vitello, P.
        (2000).
        JWL++: A Simple Reactive Flow Code Package for Detonation
        Propellants, Explosives, Pyrotechnics, 25, 54-58.
    \endverbatim

SourceFiles
    Murnaghan.C

\*---------------------------------------------------------------------------*/

#ifndef MurnaghanBlast_H
#define MurnaghanBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Murnaghan Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Murnaghan
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- Reference pressure
    scalar pRef_;

    //- Coefficients
    scalar n_;
    scalar kappa_;

    //- Mie Grunsien coefficient
    scalar Gamma_;


public:

    // Constructors

        //- Multiphase constuctor
        Murnaghan(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "Murnaghan<"
              + word(Specie::typeName()) + '>';
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return true;
        }

        //- Return pressure
        scalar p
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Gamma function
        scalar Gamma
        (
            const scalar& rho,
            const scalar& T,
            const scalar& cv
        ) const;

        //- Speed of sound
        scalar speedOfSound
        (
            const scalar& p,
            const scalar& rho,
            const scalar& T,
            const scalar& cv
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar& p,
            const scalar& rho,
            const scalar& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MurnaghanI.H"

#ifdef NoRepository
    #include "Murnaghan.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
