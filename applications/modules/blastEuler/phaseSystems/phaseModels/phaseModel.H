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
    Foam::phaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    phaseModel.C
    phaseModelNew.C
    phaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "diameterModel.H"
#include "fvc.H"
#include "fvm.H"
#include "timeIntegrationSystem.H"
#include "runTimeSelectionTables.H"
#include "viscosity.H"
#include "compressibleMomentumTransportModel.H"
#include "phaseFluidThermophysicalTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "fluidBlastThermo.H"
#include "solidBlastThermo.H"
#include "hashedWordList.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class phaseSystem;
class phaseFluxScheme;

// Trait for converting the ThermoModel's thermo type to the thermo type needed
// for the thermophysical transport model type; i.e., from rho-type thermo to
// fluid-type thermo.

template<class ThermoModel>
struct PhaseModelTransportThermoModel;

template<>
struct PhaseModelTransportThermoModel<fluidBlastThermo>
{
    typedef fluidThermo type;
};


/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/


class phaseModel
:
    public volScalarField,
    public timeIntegrationSystem,
    public viscosity
{
protected:
    // Protected typedefs

        //- Thermo type for the thermophysical transport model
        typedef
            typename PhaseModelTransportThermoModel
            <
                fluidBlastThermo
            >::type
            transportThermoModel;


    // Protected data

        //- Reference to phase system
        const phaseSystem& fluid_;

        //- Name of phase
        word name_;

        //- Index of phase
        label index_;

        //- Phase dictionary
        dictionary phaseDict_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;


        //- Primitative variables

            //- Mean velocity
            volVectorField U_;


        // Conseritative variables

            //- Phasic density
            volScalarField alphaRho_;

            //- Phasic momentum
            volVectorField alphaRhoU_;

            //- Phasic total energy
            volScalarField alphaRhoE_;


        // Conservative fluxes

            //- Volumetric flux of the phase
            surfaceScalarField phi_;

            //- Phasic mass flux
            surfaceScalarField alphaRhoPhi_;

            //- Phasic momentum flux
            surfaceVectorField alphaRhoUPhi_;

            //- Phasic total energy flux
            surfaceScalarField alphaRhoEPhi_;


        //- Diameter model
        UautoPtr<diameterModel> dPtr_;

        //- Turbulence model
        autoPtr<phaseCompressible::momentumTransportModel> turbulence_;

        //- Thermophysical transport model
        autoPtr
        <
            PhaseThermophysicalTransportModel
            <
                phaseCompressible::momentumTransportModel,
                transportThermoModel
            >
        > thermophysicalTransport_;

        //- Is the volume fraction transported
        bool solveAlpha_;

        //- Solution vector
        vector solutionDs_;

        // Protected functions
        virtual tmp<volScalarField> ESource() const = 0;

        //- Solve phase density
        //  Seperated to allow solving before thermodynamics,
        //  and solve thermo before energy
        void solveAlphaRho();

        //- Solve the diameter model
        virtual void solveD();


public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            dictionary,
            (
                const phaseSystem& fluid,
                const word& phaseName,
                const label index
            ),
            (fluid, phaseName, index)
        );

    // Constructors
        phaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const phaseSystem& fluid_;
            mutable label indexCounter_;

        public:

            iNew
            (
                const phaseSystem& fluid
            )
            :
                fluid_(fluid),
                indexCounter_(-1)
            {}

            autoPtr<phaseModel> operator()(Istream& is) const
            {
                indexCounter_++;
                return autoPtr<phaseModel>
                (
                    phaseModel::New(fluid_, word(is), indexCounter_)
                );
            }
        };

         //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        using timeIntegrationSystem::mesh;

        //- Return that the phase uses a slave pressure
        virtual bool slavePressure() const
        {
            return false;
        }

        //- Does the phase model track total energy
        virtual bool totalEnergy() const
        {
            return true;
        }

        //- Is the phase granular
        virtual bool granular() const
        {
            return false;
        }

        //- Return phase system
        const phaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Determine if alpha is solved
        virtual void solveAlpha(const bool s);

        //- Determine if alpha is solved
        bool solveAlpha() const
        {
            return solveAlpha_;
        }

        //- Return the name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the name of this phase
        const word& keyword() const
        {
            return name_;
        }

        //- Return the phase dictinoary
        const dictionary& dict() const
        {
            return phaseDict_;
        }

        //- Return phase index
        label index() const
        {
            return index_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return the solution directions
        const vector& solutionDs() const
        {
            return solutionDs_;
        }

        //- Return the flux scheme
        virtual const phaseFluxScheme& flux() const = 0;

        //- Access the flux scheme
        virtual phaseFluxScheme& flux() = 0;

        //- Turbulence model
        const phaseCompressible::momentumTransportModel& turbulence() const
        {
            return turbulence_();
        }

        //- Thermophysical transport model
        const
            PhaseThermophysicalTransportModel
            <
                phaseCompressible::momentumTransportModel,
                transportThermoModel
            >&
        thermophysicalTransport() const
        {
            return thermophysicalTransport_();
        }

        //- Return the residual phase-fraction for given phase
        virtual dimensionedScalar residualAlpha() const
        {
            return thermo().residualAlpha();
        }

        //- Return the residual phase-dsnity for given phase
        virtual dimensionedScalar residualRho() const
        {
            return thermo().residualRho();
        }

        //- Return the residual phase-density for given phase
        virtual dimensionedScalar residualAlphaRho() const
        {
            return residualAlpha()*residualRho();
        }


        // integrationSystem functions

            //- Scale volume fraction
            virtual void scaleVolumeFraction
            (
                const scalar sumAlpha,
                const label celli
            );

            //- Correct volume fraction
            virtual void correctVolumeFraction
            (
                const scalar alpha,
                const label celli
            );

            //- Decode conservative variables
            virtual void decode() = 0;

            //- Encode conservative variables
            virtual void encode() = 0;

            //- Update fluxes
            virtual void update();

            //- Solve sub-step stepi
            virtual void solve() = 0;

            //- Solve viscous and source terms
            virtual void postUpdate() = 0;

            //- Correct delta terms
            virtual void correctDeltas()
            {}


        //- Return conserved quantities

            //- Return the density
            virtual const volScalarField& rho() const = 0;

            //- Access the density
            virtual volScalarField& rho() = 0;

            //- Constant access to the mean velocity
            virtual const volVectorField& U() const
            {
                return U_;
            }

            //- Non-const access to the mean velocity
            virtual volVectorField& U()
            {
                return U_;
            }

            //- Constant access to granular temperature
            virtual tmp<volScalarField> Theta() const
            {
                NotImplemented;
                return alphaRho_;
            }

            //- Non-constant access to granular temperature
            virtual volScalarField& Theta()
            {
                NotImplemented;
                return alphaRho_;
            }

            //- Return diameterModel
            virtual const diameterModel& dModel() const
            {
                return dPtr_();
            }

            //- Return diameter of nodei
            virtual tmp<volScalarField> d() const
            {
                return dPtr_->d();
            }

            //- Return diameter of nodei for celli
            virtual scalar celld(const label celli) const
            {
                return dPtr_->d()[celli];
            }

            //- Return the thermophysical pressure
            virtual const volScalarField& p() const = 0;

            //- Return non-const access to the thermodynamic pressure
            virtual volScalarField& p() = 0;

            //- Return non-const access to the pressure gradient
            virtual tmp<volVectorField> gradP() const = 0;

            //- Return non-const access to the volume fraction gradient
            virtual tmp<volVectorField> gradAlpha() const = 0;

            //- Return the internal energy
            virtual const volScalarField& he() const
            {
                return thermo().he();
            }

            //- Return non-const access to the internal energy
            virtual volScalarField& he()
            {
                return thermo().he();
            }

            //- Return temperature
            virtual const volScalarField& T() const
            {
                return thermo().T();
            }

            //- Return surface temperature
            virtual tmp<volScalarField> Ts() const
            {
                return thermo().T();
            }

            //- Constant access to phasic density
            virtual const volScalarField& alphaRho() const
            {
                return alphaRho_;
            }

            //- Non-constant access to phasic density
            virtual volScalarField& alphaRho()
            {
                return alphaRho_;
            }

            //- Constant access to phasic momentum
            virtual const volVectorField& alphaRhoU() const
            {
                return alphaRhoU_;
            }

            //- Non-constant access to phasic momentum
            virtual volVectorField& alphaRhoU()
            {
                return alphaRhoU_;
            }

            //- Constant access to phasic total energy
            virtual const volScalarField& alphaRhoE() const
            {
                return alphaRhoE_;
            }

            //- Non-constant access to phasic total energy
            virtual volScalarField& alphaRhoE()
            {
                return alphaRhoE_;
            }

            //- Constant access to phasic pseudo thermal energy
            virtual tmp<volScalarField> alphaRhoPTE() const
            {
                NotImplemented;
                return alphaRhoE_;
            }

            //- Non-constant access to phasic pseudoo thermal energy
            virtual volScalarField& alphaRhoPTE()
            {
                NotImplemented
                return alphaRhoE_;
            }


        // Return conservative fluxes

            //- Constant access to the volumetric flux
            tmp<surfaceScalarField> phi() const
            {
                return phi_;
            }

            //- Constant access to the volume fraction flux
            virtual tmp<surfaceScalarField> alphaPhi() const;

            //- Constant access to the mass flux
            virtual tmp<surfaceScalarField> alphaRhoPhi() const
            {
                return alphaRhoPhi_;
            }


        // Thermodynamic properties

            //- Return thermodynamic model
            virtual const blastThermo& thermo() const = 0;

            //- Access thermodynamic model
            virtual blastThermo& thermo() = 0;

            //- Return thermodynamic model
            virtual const ::Foam::fluidThermo& fluidThermo() const = 0;

            //- Return the laminar viscosity
            virtual tmp<volScalarField> nu() const = 0;

            //- Return the laminar viscosity for patch
            virtual tmp<scalarField> nu(const label patchi) const = 0;

            //- Return the laminar viscosity for patch
            virtual scalar cellnu(const label celli) const = 0;

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const = 0;

            //- Return the laminar dynamic viscosity for patch
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Return the speed of sound
            virtual tmp<volScalarField> speedOfSound() const = 0;

            //- Return cell pressure
            virtual scalar cellpRhoT(const label celli) const = 0;

            //- Return cell speedOfSound
            virtual scalar cellSpeedOfSound
            (
                const scalar p,
                const label celli
            ) const = 0;


        //- Return the PTE production
        virtual tmp<volScalarField> productionSource(const phaseModel&) const
        {
            NotImplemented;
            return *this;
        }

        //- Return the PTE dissipation rate
        virtual tmp<volScalarField>
        dissipationSource(const phaseModel&, const dimensionedScalar&) const
        {
            NotImplemented;
            return *this;
        }

        //- Correct
        virtual void correct()
        {}

        //- Solve change in size moments based on breakup and coalesence
        virtual void solveSource()
        {
            return;
        }

        //- Correct thermodynamic models
        virtual void correctThermo() = 0;

        //- Return the maximum realizable courant number
        virtual scalar realizableCo() const
        {
            return 1.0;
        }

        //- Return the max courant number
        virtual scalar CoNum() const
        {
            return 0.0;
        }

        //- Read phase properties dictionary
        virtual bool read()
        {
            return true;
        }

        //- Read phase properties dictionary
        virtual bool read(const bool readOK)
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
