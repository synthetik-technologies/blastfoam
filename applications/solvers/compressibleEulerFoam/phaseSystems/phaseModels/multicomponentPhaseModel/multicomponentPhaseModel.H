/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::multicomponentPhaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    multicomponentPhaseModel.C
    newmulticomponentPhaseModel.C
    multicomponentPhaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentPhaseModel_H
#define multicomponentPhaseModel_H

#include "phaseModel.H"
#include "fluidThermoModel.H"
#include "fluxScheme.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class multicomponentPhaseModel Declaration
\*---------------------------------------------------------------------------*/


class multicomponentPhaseModel
:
    public phaseModel
{
protected:
    // Protected data

        //- Names of components
        wordList components_;

        //- Phase pressure
        volScalarField p_;

        //- Density field
        volScalarField rho_;

        //- Volume fraction flux
        surfaceScalarField alphaPhi_;

        //- Equation of state for phases
        PtrList<fluidThermoModel> thermos_;

        //- Phase Volume fraction of components
        PtrList<volScalarField> Ys_;

        //- Residual alpha
        dimensionedScalar residualAlpha_;

        //- Volume weighted density of phases
        PtrList<volScalarField> alphaRhoYs_;

        //- Mass fluxes
        PtrList<surfaceScalarField> alphaRhoYPhis_;

        //- Old values for ode solver
        PtrList<volScalarField> alphaOld_;
        PtrList<PtrList<volScalarField>> alphaRhoYsOld_;

        //- Stored delta fields
        PtrList<volScalarField> deltaAlpha_;
        PtrList<PtrList<volScalarField>> deltaAlphaRhoYs_;

        //- Flux scheme
        autoPtr<fluxScheme> fluxScheme_;

        //- Turbulence model
        autoPtr<PhaseCompressibleTurbulenceModel<phaseModel>> turbulence_;

        virtual tmp<volScalarField> ESource() const;


public:

    //- Runtime type information
    ClassName("multicomponent");


    // Constructor
    multicomponentPhaseModel
    (
        const phaseSystem& fluid,
        const word& phaseName,
        const label index
    );


    //- Destructor
    virtual ~multicomponentPhaseModel();


    // Member Functions

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha() const
        {
            return residualAlpha_;
        }


        // integrationSystem functions

             virtual void solveAlpha(const bool s);

            //- Decode primative variables
            virtual void decode();

            //- Encode conserved variables
            virtual void encode();

            //- Update fluxes
            virtual void update();

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


        //- Return the density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the thermophysical pressure
        virtual const volScalarField& p() const
        {
            return p_;
        }

        //- Return non-const access to the thermodynamic pressure
        virtual volScalarField& p()
        {
            return p_;
        }

        //- Return non-const access to the pressure gradient
        virtual tmp<volVectorField> gradP() const;

        //- Return non-const access to the volume fraction gradient
        virtual tmp<volVectorField> gradAlpha() const;

        //- Constant access to the volume fraction flux
        virtual tmp<surfaceScalarField> alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Return the laminar viscosity for patch
        virtual scalar nui(const label celli) const;

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const;

        //- Return the laminar dynamic viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volScalarField> alpha() const;

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const;

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature of mixture
        //  for patch [J/m/s/K]
        virtual tmp<scalarField> kappa(const label patchi) const;

        //- Thermal diffusivity for temperature of mixture
        //  for cell [J/m/s/K]
        virtual scalar kappai(const label celli) const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField& alphat
        ) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Return the specific heat capacity at constant pressure
        virtual tmp<volScalarField> Cp() const;

        //- Return the specific heat capacity at constant volume
        virtual tmp<volScalarField> Cv() const;

        //- Return the specific heat capacity at constant volume
        virtual scalar Cvi(const label celli) const;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const;

        //- Correct thermodynamic quantities
        virtual void correctThermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
