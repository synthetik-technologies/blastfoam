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
    Foam::phaseCompressibleSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    phaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseCompressibleSystem_H
#define singlePhaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseCompressibleSystem.H"
#include "fluidThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class singlePhaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseCompressibleSystem
:
    public phaseCompressibleSystem
{
protected:
// Protected data

    //- Temperature
    volScalarField T_;

    //- Internal energy
    volScalarField e_;

    //- Equation of state
    autoPtr<fluidThermoModel> thermo_;



    // ODE variables

        //- Old values for ode solver
        PtrList<volScalarField> rhoOld_;

        //- Stored delta fields
        PtrList<volScalarField> deltaRho_;

        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho()
        {
            return;
        }

public:

    TypeName("singlePhaseCompressibleSystem");

    // Constructor
    singlePhaseCompressibleSystem(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~singlePhaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const;

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const boolList& storeFields,
            const boolList& storeDeltas
        );

        //- Remove stored fields
        virtual void clearODEFields();


    // Member Access Functions

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

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const;


    //- Return functions for transport model

        //- Return calculated temperature field
        virtual tmp<volScalarField> T() const
        {
            return T_;
        }

        //- Heat capacit at contant volume of mixture [J/kg/K]
        virtual tmp<volScalarField> Cv() const;

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

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
