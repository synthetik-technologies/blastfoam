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
    Foam::thermoModel

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    thermoI.H
    thermo.C

\*---------------------------------------------------------------------------*/

#ifndef thermoModelBlast_H
#define thermoModelBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "thermoModel.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ThermoType> class thermoModel;

template<class ThermoType>
inline thermoModel<ThermoType> operator+
(
    const thermoModel<ThermoType>&,
    const thermoModel<ThermoType>&
);

template<class ThermoType>
inline thermoModel<ThermoType> operator*
(
    const scalar,
    const thermoModel<ThermoType>&
);

template<class ThermoType>
inline thermoModel<ThermoType> operator==
(
    const thermoModel<ThermoType>&,
    const thermoModel<ThermoType>&
);


/*---------------------------------------------------------------------------*\
                           Class thermoModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class thermoModel
:
    public ThermoType
{
protected:

    // Protected data

        //- Convergence criteria for temperature
        scalar tolerance_;

        //- Maximum number of iterations
        label maxIter_;


        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(rho, T)
        //  and dF(rho, T)/dT
        inline scalar T
        (
            const scalar f,
            const scalar e,
            const scalar rho,
            const scalar T0,
            scalar (thermoModel::*F)(const scalar, const scalar, const scalar) const,
            scalar (thermoModel::*dFdT)(const scalar, const scalar, const scalar) const,
            scalar (thermoModel::*limit)(const scalar) const
        ) const;


public:

    //- The thermodynamics of the individual species'
    typedef thermoModel<ThermoType> thermoType;

    // Constructors

        //- Construct from components
        inline thermoModel(const ThermoType& sp);

        //- Construct from dictionary
        thermoModel(const dictionary& dict);

        //- Construct as named copy
        inline thermoModel(const word& name, const ThermoType&);


    //- Return the instantiated type name
        static word typeName()
        {
            return  word(ThermoType::typeName());
        }


    // Member Functions

        //- Return Mie Gruniesen coefficient
        scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return pressure
        scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdRho
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return derivative of pressure w.r.t internal energy
        scalar dpde
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return speed of sound
        scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return specific heat ratio
        scalar CpByCv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Initialize internal energy
        scalar initializeEnergy
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Initialize internal energy
        scalar initializeRho
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;


    // Mass specific derived properties

            //- Gibbs free energy [J/kg]
            inline scalar G
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;

            //- Helmholtz free energy [J/kg]
            inline scalar A
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;


    // Equilibrium reaction thermodynamics

            //- Equilibrium constant [] i.t.o fugacities
            //  = PIi(fi/Pstd)^nui
            inline scalar K
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar p,
                const scalar rho,
                const scalar T,
                const scalar n
            ) const;

        // Mole specific derived properties

            //- Heat capacity at constant pressure [J/kmol/K]
            inline scalar cp
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/kmol/K]
            inline scalar s
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Enthalpy/Internal energy [J/kmol]
            inline scalar he
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Heat capacity at constant volume [J/kmol/K]
            inline scalar cv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Sensible internal energy [J/kmol]
            inline scalar es
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Absolute internal energy [J/kmol]
            inline scalar ea
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Gibbs free energy [J/kmol]
            inline scalar g
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;

            //- Helmholtz free energy [J/kmol]
            inline scalar a
            (
                const scalar p,
                const scalar rho,
                const scalar T
            ) const;


        // Derivative term used for Jacobian

            //- Derivative of B (acooding to Niemeyer et al.) w.r.t. temperature
            inline scalar dKcdTbyKc(const scalar p, const scalar T) const;

            //- Derivative of cp w.r.t. temperature
            inline scalar dcpdT(const scalar p, const scalar T) const;


    // Energy->temperature  inversion functions

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar TRhoE
            (
                const scalar T,
                const scalar rho,
                const scalar E
            ) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa
            (
                const scalar Ha,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs
            (
                const scalar Es,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa
            (
                const scalar Ea,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const;


    // Member Operators

        inline void operator+=(const thermoModel&);
        inline void operator*=(const scalar);


    // Friend operators

        friend thermoModel operator+ <ThermoType>
        (
            const thermoModel&,
            const thermoModel&
        );

        friend thermoModel operator* <ThermoType>
        (
            const scalar s,
            const thermoModel&
        );

        friend thermoModel operator== <ThermoType>
        (
            const thermoModel&,
            const thermoModel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModelI.H"

#ifdef NoRepository
    #include "thermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineThermoType(Transport, Thermo, EoS, Specie)                    \
typedef                                                                     \
    Transport                                                               \
    <                                                                       \
        thermoModel                                                         \
        <                                                                   \
            Thermo                                                          \
            <                                                               \
                EoS                                                         \
                <                                                           \
                    Specie                                                  \
                >                                                           \
            >                                                               \
        >                                                                   \
    > Transport##Thermo##EoS##Specie;

#define defineThermoTypes(Transport, EoS, Specie)                           \
    defineThermoType(Transport, eConst, EoS, Specie)                        \
    defineThermoType(Transport, hConst, EoS, Specie)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
