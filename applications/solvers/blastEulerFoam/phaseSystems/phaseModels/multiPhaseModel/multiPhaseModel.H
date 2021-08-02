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
    Foam::multiPhaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    multiPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiPhaseModel_H
#define multiPhaseModel_H

#include "phaseModel.H"
#include "fluxScheme.H"
#include "multiphaseFluidBlastThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class multiPhaseModel Declaration
\*---------------------------------------------------------------------------*/


class multiPhaseModel
:
    public phaseModel
{
protected:
    // Protected data

        //- Names of components
        wordList phases_;

        //- Phase pressure
        volScalarField p_;

        //- Density field
        volScalarField rho_;

        //- Volume fraction flux
        surfaceScalarField alphaPhi_;

        //- Primitative variables

            //- Equation of states
            multiphaseFluidBlastThermo thermo_;

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


        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho();

        //- Flux scheme
        autoPtr<fluxScheme> fluxScheme_;

        virtual tmp<volScalarField> ESource() const;


public:

    //- Runtime type information
    ClassName("multiphase");


    // Constructor
    multiPhaseModel
    (
        const phaseSystem& fluid,
        const word& phaseName,
        const label index
    );


    //- Destructor
    virtual ~multiPhaseModel();


    // Member Functions

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha() const
        {
            return thermo_.residualAlpha();
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
            virtual void solve();

            //- Update viscous terms
            virtual void postUpdate();

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

        //- Return thermodynamic model
        virtual const basicBlastThermo& thermo() const
        {
            return thermo_;
        }

        //- Access thermodynamic model
        virtual basicBlastThermo& thermo()
        {
            return thermo_;
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
