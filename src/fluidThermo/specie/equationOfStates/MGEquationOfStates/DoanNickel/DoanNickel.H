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
    Foam::equationOfStates::DoanNickel

Description
    Real equation of state for a gas

    References:
    \verbatim
        Needham, C.E., Troutman, W.W. (1966).
        THEORETICAL CALCULATIONS OF THE PHENOMENOLOGY OF HE
        DETONATIONS. VOLUME 1
        Air Force Weapons Laboratory
    \endverbatim

SourceFiles
    DoanNickelI.H
    DoanNickel.C

\*---------------------------------------------------------------------------*/

#ifndef DoanNickel_H
#define DoanNickel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DoanNickel Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class DoanNickel
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_ = 1.293;

    //- Gamma coefficients
    scalar a_ = 0.161;
    scalar b_ = 0.255;
    scalar c_ = 0.280;
    scalar d_ = 0.137;
    scalar g_ = 0.050;

    scalar a1_ = 0.048;
    scalar a2_ = 0.032;
    scalar a3_ = 0.045;

    //- Equation of state coefficients for first regime
    //  O_2 -> 2O
    scalar E1_ = 4.46;
    scalar E11_ = 0.357;
    scalar E1Offset_ = 8.5;
    scalar deltaE1_ = 0.975;
    scalar deltaE1Pow_ = 0.05;

    //- Equation of state coefficients for second regime
    //  N_2 -> 2N
    scalar E2_ = 6.63;
    scalar E22_ = 45.0;
    scalar deltaE2_ = 4.0;
    scalar E2Pow_ = 0.0157;
    scalar deltaE2Pow_ = 0.085;

    //- Equation of state coefficients for third regime
    //  electrons
    scalar E3_ = 25.5;
    scalar E33_ = 160.0;
    scalar deltaE3_ = 6.0;

    //- Max value used in exp function
    const scalar expMax_ = 250;

    //- Derivative of specific heat ration w.r.t. density
    scalar dGammadRho
    (
        const scalar& rho,
        const scalar& e
    ) const;

    //- Derivative of specific heat w.r.t. internal energy
    scalar dGammade
    (
        const scalar& rho,
        const scalar& e
    ) const;


public:

    // Constructors

        //- Single phase constructor
        DoanNickel(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "DoanNickel<"
              + word(Specie::typeName()) + '>';
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return false;
        }

        //- Return pressure deviation from ideal
        scalar Pi
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Return speed of sound
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
            const scalar& e
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpde
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& e
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

#include "DoanNickelI.H"

#ifdef NoRepository
    #include "DoanNickel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

