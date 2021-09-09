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
    Foam::phaseCompressibleSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    phaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef psiuCompressibleSystem_H
#define psiuCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "integrationSystem.H"
#include "fluxScheme.H"
#include "psiuReactionThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "uniformDimensionedFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class psiuCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class psiuCompressibleSystem
:
    public integrationSystem
{
protected:
// Protected data

    //- Thermo
    autoPtr<psiuReactionThermo> thermo_;

    //- Primitative variables

        //- Total mass
        volScalarField rho_;

        //- Velocity
        volVectorField U_;

        //- Pressure
        volScalarField& p_;

        //- Temperature
        volScalarField& T_;

        //- Internal energy
        volScalarField& e_;

        //- Refernce to unburnt internal energy
        volScalarField& eu_;

    //- Conserved variables

        //- Momentum
        volVectorField rhoU_;

        //- Energy
        volScalarField rhoE_;

        //- Total unbunrt energy
        volScalarField rhoEu_;


    //- Fluxes

        //- Volumetric flux
        surfaceScalarField phi_;

        //- Total mass flux
        surfaceScalarField rhoPhi_;

        //- Momentum flux
        surfaceVectorField rhoUPhi_;

        //- Energy flux
        surfaceScalarField rhoEPhi_;

        //- Flux scheme
        autoPtr<fluxScheme> fluxScheme_;

        //- Turbulence model
        autoPtr<compressible::momentumTransportModel> turbulence_;

        //- Thermo turbulence model
        autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport_;


        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;


public:

    TypeName("psiuCompressibleSystem");

    // Constructor
    psiuCompressibleSystem(const fvMesh& mesh);

    //- Destructor
    virtual ~psiuCompressibleSystem();


    // Member functions

        //- Non-const reference to thermo
        psiuReactionThermo& thermo()
        {
            return thermo_();
        }

        //- Const reference to thermo
        const psiuReactionThermo& thermo() const
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
        virtual void solve();

        //- Set old lists and fluxes (initialization of fields)
        virtual void postUpdate();


    // Member Access Functions

        //- Return the flux scheme
        virtual const fluxScheme& flux() const
        {
            return fluxScheme_();
        }

        //- Access the flux scheme
        virtual fluxScheme& flux()
        {
            return fluxScheme_();
        }

        //- Const access to turbulence model
        const compressible::momentumTransportModel& turbulence() const
        {
            return turbulence_();
        }

        //- Non-const access to turbulence model
        compressible::momentumTransportModel& turbulence()
        {
            return turbulence_();
        }

        //- Const access to fluidThermophysicalTransportModel model
        const fluidThermophysicalTransportModel& thermophysicalTransport() const
        {
            return thermophysicalTransport_();
        }

        //- Non-const access to fluidThermophysicalTransportModel model
        fluidThermophysicalTransportModel& thermophysicalTransport()
        {
            return thermophysicalTransport_();
        }

        //- Velocity
        virtual tmp<volVectorField> U() const;

        //- Const access to internal energy
        virtual const volScalarField& e() const
        {
            return e_;
        }

        //- Non-const access to internal energy
        virtual volScalarField& e()
        {
            return e_;
        }

        //- Return volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return mass flux
        const surfaceScalarField& rhoPhi() const
        {
            return rhoPhi_;
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

        //- Density of mixture [kg/m^3]
        virtual tmp<volScalarField> rho() const;

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
