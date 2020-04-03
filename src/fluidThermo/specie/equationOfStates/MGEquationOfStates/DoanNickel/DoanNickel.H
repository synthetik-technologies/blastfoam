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

    Example entry for air from the Needham and Troutman (1966):
    \verbatim
    mixture
    {
        type        DoanNickel;

        rho0 1.293;
        e0 206730.6122;

        a 0.161;
        b 0.256;
        c 0.280;
        d 0.137;
        g 0.5;

        a1 0.161;
        a2 0.32;
        a3 0.015;

        E1 4.46e6;
        E2 6.63e6;
        E3 25.5e6;

        E11 0.357e6;
        E1Offset 8.5e6;
        deltaE1 0.975e6;
        deltaE1Pow 0.05;

        E22 45.0e6;
        E2Pow 0.015;
        deltaE2 4.0e6;
        deltaE2Pow 0.095;

        E33 160.0e6;
        deltaE3 6.0e6;
    }
    \endverbatim

SourceFiles
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
    scalar rho0_;

    //- Reference internal energy
    scalar e0_;

    //- Gamma coefficients
    scalar a_;
    scalar b_;
    scalar c_;
    scalar d_;
    scalar g_;

    scalar a1_;
    scalar a2_;
    scalar a3_;

    //- Equation of state coefficients for first regime
    scalar E1_;
    scalar E11_;
    scalar E1Offset_;
    scalar deltaE1_;
    scalar deltaE1Pow_;

    //- Equation of state coefficients for second regime
    scalar E2_;
    scalar E22_;
    scalar deltaE2_;
    scalar E2Pow_;
    scalar deltaE2Pow_;

    //- Equation of state coefficients for third regime
    scalar E3_;
    scalar E33_;
    scalar deltaE3_;

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

