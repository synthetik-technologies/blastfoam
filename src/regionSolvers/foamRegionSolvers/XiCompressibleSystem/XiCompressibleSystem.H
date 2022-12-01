/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::XiCompressibleSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    XiCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef XiCompressibleSystem_H
#define XiCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystem.H"
#include "fluxScheme.H"
#include "psiuReactionThermo.H"
#include "uniformDimensionedFields.H"
#include "laminarFlameSpeed.H"
#include "ignition.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class XiCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class XiCompressibleSystem
:
    public compressibleSystem
{
protected:
// Protected data

    //- Thermo
    autoPtr<psiuReactionThermo> thermo_;

    //- Primitative variables

        //- Total mass
        volScalarField rho_;

        //- Pressure
        volScalarField& p_;

        //- Temperature
        volScalarField& T_;

        //- Internal energy
        volScalarField& e_;

        //- Refernce to unburnt internal energy
        volScalarField& eu_;


    // Additional conserved variables

        //- Total unbunrt energy
        volScalarField rhoEu_;


        //- Composition
        basicCombustionMixture& composition_;

        //- Reaction progress variable
        volScalarField& b_;


    // Flame models

        //-
        volScalarField Xi_;

        // Flame speed model
        autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed_;

        // Laminar flame speed
        volScalarField Su_;

        //- Turbulent flame speed
        volScalarField St_

        //- Bounds
        dimensionedScalar SuMin_;
        dimensionedScalar SuMax_;


public:

    TypeName("Xi");

    // Constructor
    XiCompressibleSystem(const fvMesh& mesh);

    //- Destructor
    virtual ~XiCompressibleSystem();


    // Member functions

        //- Non-const reference to reaction thermo
        psiuReactionThermo& rThermo()
        {
            return thermo_();
        }

        //- Const reference to reaction thermo
        const psiuReactionThermo& rThermo() const
        {
            return thermo_();
        }

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual bool solve();

        //- Set old lists and fluxes (initialization of fields)
        virtual void postUpdate();


    // Member Access Functions

        //- Const access to thermo model
        virtual const fluidThermo& thermo() const
        {
            return thermo_();
        }

        //- Non-const access to thermo model
        virtual fluidThermo& thermo()
        {
            return thermo_();
        }

        //- Density of mixture [kg/m^3]
        virtual const volScalarField& rho() const;

        //- Const access to internal energy
        virtual const volScalarField& he() const
        {
            return e_;
        }

        //- Non-const access to internal energy
        virtual volScalarField& he()
        {
            return e_;
        }

        //- Return thermodynamic pressure
        const volScalarField& p() const
        {
            return p_;
        }

        //- Return temperature
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const;


    //- Return functions for transport model

        //- Density of unburnt mixture [kg/m^3]
        virtual tmp<volScalarField> rhou() const;

        //- Heat capacit at contant volume of mixture [J/kg/K]
        virtual tmp<volScalarField> Cv() const;

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Dynamic viscosity of unburnt mixture [kg/m/s]
        virtual tmp<volScalarField> muu() const;

        //- Dynamic viscosity of burnt mixture [kg/m/s]
        virtual tmp<volScalarField> mub() const;

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volScalarField> alpha() const;

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff(const volScalarField&) const;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const;

        //- Thermal diffusivity for temperature of mixture [W/m/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappa(const label patchi) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff(const volScalarField&) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
