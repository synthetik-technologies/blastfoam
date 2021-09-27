/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied technologies
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
    Foam::simpleBlastThermo

Description
    A simple thermo class used to define a runTimeSelection table for testing
    thermodynamic models

SourceFiles
    simpleBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef simpleBlastThermo_H
#define simpleBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "wordIOList.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class simpleBlastThermo
{
//- Private functions

    //- Generic lookup for thermodynamics package thermoTypeName
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoTypeDict,
            Table* tablePtr,
            const int nCmpt,
            const char* cmptNames[],
            const word& thermoTypeName
        );

        //- Generic lookup for each of the related thermodynamics packages
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoDict,
            Table* tablePtr
        );

        //- Return the thermo type name given a dictionary
        static word readThermoType(const dictionary& dict);

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName
        );

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName,
            const label nCmpts
        );
public:

    //- Runtime type information
    TypeName("simpleBlastThermo");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            simpleBlastThermo,
            dictionary,
            (const dictionary& dict),
            (dict)
        );

    // Constructor

        //- Construct with input field
        simpleBlastThermo();


    // Selectors

        static autoPtr<simpleBlastThermo> New(const dictionary& dict);

    //- Destructor
    virtual ~simpleBlastThermo();

    // Access functions

        //- Species

            virtual scalar W() const = 0;
            virtual scalar Y() const = 0;
            virtual scalar R() const = 0;

        //- Equation of state

            //- Return derivative of pressure w.r.t. specific volume
            virtual scalar dpdv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return derivative of pressure w.r.t. Temperature
            virtual scalar dpdT
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return derivative of pressure w.r.t. internal energy
            virtual scalar dpde
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Internal energy correction [J/kg]
            virtual scalar E
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Enthalpy correction [J/kg]
            virtual scalar H
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat capacity difference [J/kg/K]
            virtual scalar CpMCv
            (
                const scalar rho,
                const scalar e,
                const scalar T,
                const scalar CpCv
            ) const = 0;

        //- Thermo model

            //- Heat capacity at constant volume of mixture [J/kg/K]
            virtual scalar Cv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Absolute internal energy [J/kg]
            virtual scalar Ea
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure of mixture [J/kg/K]
            virtual scalar Cp
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat of formation [J/kg]
            virtual scalar Hf() const = 0;

            //- Flame temperature [K]
            virtual scalar flameT() const = 0;


            //- Entropy [J/kg/K]
            virtual scalar S
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Derivative of Gibbs free energy w.r.t. temperature
            virtual scalar dGdT(const scalar p, const scalar T) const = 0;

            //- Temperature derivative of heat capacity at constant pressure
            virtual scalar dCpdT(const scalar p, const scalar T) const = 0;


        //- Master thermo class

            //- Return Mie Gruniesen coefficient
            virtual scalar Gamma
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return pressure
            virtual scalar p
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return speed of sound
            virtual scalar cSqr
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return specific heat ratio
            virtual scalar CpByCv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Initialize internal energy
            virtual scalar initializeEnergy
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Initialize internal energy
            virtual scalar rhoPT
            (
                const scalar p,
                const scalar T
            ) const = 0;

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            virtual scalar TRhoE
            (
                const scalar T,
                const scalar rho,
                const scalar E
            ) const = 0;

            //- Temperature from sensible enthalpy given an initial T0
            virtual scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const = 0;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            virtual scalar THa
            (
                const scalar Ha,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const = 0;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            virtual scalar TEs
            (
                const scalar Es,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const = 0;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            virtual scalar TEa
            (
                const scalar Ea,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const = 0;


        // Transport

            //- Dynamic viscosity [kg/ms]
            virtual scalar mu
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/mK]
            virtual scalar kappa
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Thermal diffusivity of enthalpy [kg/ms]
            virtual scalar alphah
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;
};


template<class ThermoType>
class SimpleBlastThermo
:
    public ThermoType,
    public simpleBlastThermo
{
public:

    //- Runtime type information
    TypeName("SimpleBlastThermo");

    //- Construct with input field
    SimpleBlastThermo(const dictionary& dict)
    :
        ThermoType(dict)
    {}

    //- Destructor
    virtual ~SimpleBlastThermo()
    {}

    //- Species

            virtual scalar W() const
            {
                return ThermoType::W();
            }

            virtual scalar Y() const
            {
                return ThermoType::Y();
            }

            virtual scalar R() const
            {
                return ThermoType::R();
            }

        //- Equation of state

            //- Return derivative of pressure w.r.t. specific volume
            virtual scalar dpdv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::dpdv(rho, e, T);
            }

            //- Return derivative of pressure w.r.t. Temperature
            virtual scalar dpdT
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::dpdT(rho, e, T);
            }

            //- Return derivative of pressure w.r.t. internal energy
            virtual scalar dpde
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::dpde(rho, e, T);
            }

            //- Internal energy correction [J/kg]
            virtual scalar E
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::E(rho, e, T);
            }

            //- Enthalpy correction [J/kg]
            virtual scalar H
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::H(rho, e, T);
            }

            //- Heat capacity difference [J/kg/K]
            virtual scalar CpMCv
            (
                const scalar rho,
                const scalar e,
                const scalar T,
                const scalar CpCv
            ) const
            {
                return ThermoType::CpMCv(rho, e, T, CpCv);
            }

        //- Thermo model

            //- Heat capacity at constant volume of mixture [J/kg/K]
            virtual scalar Cv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Cv(rho, e, T);
            }

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Es(rho, e, T);
            }

            //- Absolute internal energy [J/kg]
            virtual scalar Ea
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Ea(rho, e, T);
            }

            //- Heat capacity at constant pressure of mixture [J/kg/K]
            virtual scalar Cp
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Cp(rho, e, T);
            }

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Hs(rho, e, T);
            }

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Ha(rho, e, T);
            }

            //- Heat of formation [J/kg]
            virtual scalar Hf() const
            {
                return ThermoType::Hf();
            }

            //- Flame temperature [K]
            virtual scalar flameT() const
            {
                return ThermoType::flameT();
            }


            //- Entropy [J/kg/K]
            virtual scalar S
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::S(p, rho, e, T);
            }

            //- Derivative of Gibbs free energy w.r.t. temperature
            virtual scalar dGdT(const scalar p, const scalar T) const
            {
                return ThermoType::dGdT(p, T);
            }

            //- Temperature derivative of heat capacity at constant pressure
            virtual scalar dCpdT(const scalar p, const scalar T) const
            {
                return ThermoType::dCpdT(p, T);
            }


        //- Master thermo class

            //- Return Mie Gruniesen coefficient
            virtual scalar Gamma
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::Gamma(rho, e, T);
            }

            //- Return pressure
            virtual scalar p
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::p(rho, e, T);
            }

            //- Return speed of sound
            virtual scalar cSqr
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::cSqr(p, rho, e, T);
            }

            //- Return specific heat ratio
            virtual scalar CpByCv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::CpByCv(rho, e, T);
            }

            //- Initialize internal energy
            virtual scalar initializeEnergy
            (
                const scalar p,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::initializeEnergy(p, rho, e, T);
            }

            //- Initialize internal energy
            virtual scalar rhoPT
            (
                const scalar p,
                const scalar T
            ) const
            {
                return ThermoType::rhoPT(p,T);
            }

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            virtual scalar TRhoE
            (
                const scalar T,
                const scalar rho,
                const scalar E
            ) const
            {
                return ThermoType::TRhoE(T, rho, E);
            }

            //- Temperature from sensible enthalpy given an initial T0
            virtual scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const
            {
                return ThermoType::THs(Hs, p, rho, T0);
            }

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            virtual scalar THa
            (
                const scalar Ha,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const
            {
                return ThermoType::THa(Ha, p, rho, T0);
            }

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            virtual scalar TEs
            (
                const scalar Es,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const
            {
                return ThermoType::TEs(Es, p, rho, T0);
            }

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            virtual scalar TEa
            (
                const scalar Ea,
                const scalar p,
                const scalar rho,
                const scalar T0
            ) const
            {
                return ThermoType::TEa(Ea, p, rho, T0);
            }


        // Transport

            //- Dynamic viscosity [kg/ms]
            virtual scalar mu
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::mu(rho, e, T);
            }

            //- Thermal conductivity [W/mK]
            virtual scalar kappa
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::kappa(rho, e, T);
            }

            //- Thermal diffusivity of enthalpy [kg/ms]
            virtual scalar alphah
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                return ThermoType::alphah(rho, e, T);
            }
};


#define addSimpleEoS(Transport, Thermo, EoS)                               \
    typedef SimpleBlastThermo<Transport##Thermo##EoS##blastSpecie>                 \
        simple##Transport##Thermo##EoS;                                    \
                                                                           \
    defineTemplateTypeNameAndDebugWithName                                 \
    (                                                                      \
        simple##Transport##Thermo##EoS,                                    \
        (Transport##Thermo##EoS##blastSpecie::typeName()).c_str(),         \
        0                                                                  \
    );                                                                     \
    addToRunTimeSelectionTable                                             \
    (                                                                      \
        simpleBlastThermo,                                                         \
        simple##Transport##Thermo##EoS,                                    \
        dictionary                                                         \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "simpleBlastThermoTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
