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
    Foam::multiphaseCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    multiphaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseCompressibleSystem_H
#define multiphaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseCompressibleSystem.H"
#include "multiphaseFluidThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiphaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseCompressibleSystem
:
    public phaseCompressibleSystem
{
protected:
// Protected data

    //- Primitative variables

        //- Temperature
        volScalarField T_;

        //- Internal energy
        volScalarField e_;

        //- Equation of states
        multiphaseFluidThermo thermo_;

        //- Volume fractions
        PtrList<volScalarField>& alphas_;

        //- Phase densities
        UPtrList<volScalarField>& rhos_;


    //- Conserved variables

        //- Phase masses
        PtrList<volScalarField> alphaRhos_;


    //- Fluxes

        //- Volume fraction fluxes
        PtrList<surfaceScalarField> alphaPhis_;

        //- Phase mass fluxes
        PtrList<surfaceScalarField> alphaRhoPhis_;


    // ODE variables

        //- Old values for ode solver
        PtrList<PtrList<volScalarField>> alphasOld_;
        PtrList<PtrList<volScalarField>> alphaRhosOld_;

        //- Stored delta fields
        PtrList<PtrList<volScalarField>> deltaAlphas_;
        PtrList<PtrList<volScalarField>> deltaAlphaRhos_;

        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho();


public:

    TypeName("multiphaseCompressibleSystem");

    // Constructor
    multiphaseCompressibleSystem(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    ~multiphaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const;

        //- Solve a sub-step
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes
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
        virtual tmp<volScalarField> speedOfSound() const
        {
            return thermo_.speedOfSound();
        }

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
