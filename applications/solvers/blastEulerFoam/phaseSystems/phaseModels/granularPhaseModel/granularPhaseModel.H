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
    Foam::granularPhaseModel

Description
    Phase model to describe a solid phase

Reference:
    \verbatim
        Houim, R.W., Oran, E.S. (2016).
        A multiphase model for compressible granular–gaseous flows: formulation and
        initial tests.
        Journal of Fluid Mechanics 789, 166–220.
    \endverbatim

    \verbatim
        Huilin, L., Gidaspow, D. (2003).
        Hydrodynamics of binary fluidization in a riser: CFD simulation using two
        granular temperatures.
        Chemical Engineering Science 58, 3777–3792.
    \endverbatim

SourceFiles
    granularPhaseModel.C
    newgranularPhaseModel.C
    granularPhaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef granularPhaseModel_H
#define granularPhaseModel_H

#include "phaseModel.H"
#include "solidBlastThermo.H"
#include "kineticTheoryModel.H"
#include "phaseFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Phase>
class PhaseCompressibleTurbulenceModel;


/*---------------------------------------------------------------------------*\
                           Class granularPhaseModel Declaration
\*---------------------------------------------------------------------------*/


class granularPhaseModel
:
    public phaseModel,
    public kineticTheoryModel
{
    // Private functions

protected:
    // Protected data

        //- Thermophysical properties
        autoPtr<solidBlastThermo> thermoPtr_;

        //- Phase density
        volScalarField& rho_;

        //- Internal energy
        volScalarField& e_;

        //- Temperature
        volScalarField& T_;

        //- Pseudo thermal energy
        volScalarField alphaRhoPTE_;

        //- Pseudo thermal energy flux
        surfaceScalarField alphaRhoPTEPhi_;

        //- Flux scheme
        autoPtr<phaseFluxScheme> fluxScheme_;


public:

    //- Runtime type information
    ClassName("granular");

    // Constructors

        granularPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

        granularPhaseModel
        (
            const word& thermoType,
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~granularPhaseModel();


    //- Return that the phase uses a slave pressure
        virtual bool slavePressure() const
        {
            return true;
        }

        //- Is the phase granular
        virtual bool granular() const
        {
            return true;
        }

        //- Determine if alpha is solved
        virtual void solveAlpha(const bool s)
        {
            solveAlpha_ = false;
        }

        //- Return the residual phase-density for given phase
        virtual dimensionedScalar residualAlphaRho() const
        {
            return residualAlpha()*Foam::min(rho_);
        }


        // integrationSystem functions

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
        virtual const phaseFluxScheme& flux() const
        {
            return fluxScheme_();
        }

        //- Access the flux scheme
        virtual phaseFluxScheme& flux()
        {
            return fluxScheme_();
        }

        //- Return conserved quantities

            //- Return the density
            virtual const volScalarField& rho() const
            {
                return rho_;
            }

            //- Constant access to granular temperature
            virtual tmp<volScalarField> Theta() const
            {
                return Theta_;
            }

            //- Non-constant access to granular temperature
            virtual volScalarField& Theta()
            {
                return Theta_;
            }

            //- Return the thermophysical pressure
            virtual const volScalarField& p() const
            {
                return fluid_.p();
            }

            //- Return non-const access to the thermodynamic pressure
            virtual volScalarField& p()
            {
                NotImplemented;
                return *this;
            }

            //- Return non-const access to the pressure gradient
            virtual tmp<volVectorField> gradP() const;

            //- Return non-const access to the volume fraction gradient
            virtual tmp<volVectorField> gradAlpha() const;

            //- Constant access to phasic pseudo thermal energy
            virtual tmp<volScalarField> alphaRhoPTE() const
            {
                return alphaRhoPTE_;
            }

            //- Non-constant access to phasic pseudoo thermal energy
            virtual volScalarField& alphaRhoPTE()
            {
                return alphaRhoPTE_;
            }


        //- Return volume fraction flux
        virtual tmp<surfaceScalarField> alphaPhi() const
        {
            return alphaRhoPhi_/fvc::interpolate(rho_);
        }

        //- Return the energy source
        virtual tmp<volScalarField> ESource() const;


        // Thermodynamic properties

            //- Return thermodynamic model
            virtual const blastThermo& thermo() const
            {
                return thermoPtr_();
            }

            //- Access thermodynamic model
            virtual blastThermo& thermo()
            {
                return thermoPtr_();
            }

            //- Return the laminar viscosity
            virtual tmp<volScalarField> nu() const
            {
                return this->nut_;
            }

            //- Return the laminar viscosity for patch
            virtual tmp<scalarField> nu(const label patchi) const
            {
                return this->nut_.boundaryField()[patchi];
            }

            //- Return the laminar viscosity for patch
            virtual scalar nui(const label celli) const
            {
                return this->nut_[celli];
            }

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const
            {
                return this->nut_*rho_;
            }

            //- Return the laminar dynamic viscosity for patch
            virtual tmp<scalarField> mu(const label patchi) const
            {
                return
                    this->nut_.boundaryField()[patchi]
                   *rho_.boundaryField()[patchi];
            }

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual tmp<volScalarField> alpha() const;

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alpha(const label patchi) const;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const
            {
                return this->kappa_;
            }

            //- Thermal diffusivity for temperature of mixture
            //  for patch [J/m/s/K]
            virtual tmp<scalarField> kappa(const label patchi) const
            {
                return this->kappa_.boundaryField()[patchi];
            }

            //- Thermal diffusivity for temperature of mixture
            //  for cell [J/m/s/K]
            virtual scalar kappai(const label celli) const
            {
                return this->kappa_[celli];
            }

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const
            {
                return thermoPtr_->kappaEff(alphat);
            }

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                return thermoPtr_->kappaEff(alphat, patchi);
            }

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const
            {
                return thermoPtr_->alphaEff(alphat);
            }

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                return thermoPtr_->alphaEff(alphat, patchi);
            }

            //- Return the specific heat capacity at constant pressure
            virtual tmp<volScalarField> Cp() const;

            //- Return the specific heat capacity at constant volume
            virtual tmp<volScalarField> Cv() const;

            //- Return the specific heat capacity at constant volume
            virtual scalar Cvi(const label celli) const
            {
                return thermoPtr_->Cvi(celli);
            }

            //- Return the speed of sound
            virtual tmp<volScalarField> speedOfSound() const;


        //- Return the PTE production
        virtual tmp<volScalarField>
        productionSource(const phaseModel&) const;

        //- Return the PTE dissipation rate
        virtual tmp<volScalarField>
        dissipationSource(const phaseModel&) const;

        //- Correct
        virtual void correct()
        {}

        //- Correct thermodynamic models
        virtual void correctThermo();

        //- Read phase properties dictionary
        virtual bool read()
        {
            return true;
        }

        //- Read phase properties dictionary
        virtual bool read(const bool readOK)
        {
            return readOK;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
