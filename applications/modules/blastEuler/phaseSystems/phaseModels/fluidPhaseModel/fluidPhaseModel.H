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
    Foam::fluidPhaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    fluidPhaseModel.C
    newfluidPhaseModel.C
    fluidPhaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef fluidPhaseModel_H
#define fluidPhaseModel_H

#include "phaseModel.H"
#include "fluidBlastThermo.H"
#include "phaseFluxScheme.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class fluidPhaseModel Declaration
\*---------------------------------------------------------------------------*/


class fluidPhaseModel
:
    public phaseModel
{
protected:
    // Protected data

        //- Thermophysical properties
        autoPtr<fluidBlastThermo> thermoPtr_;

        //- Density
        volScalarField& rho_;

        //- Internal energy
        volScalarField& e_;

        //- Temperature
        volScalarField& T_;

        //- Pressure
        volScalarField& p_;

        //- Kinetic energy
        volScalarField K_;

        //- Flux scheme
        autoPtr<phaseFluxScheme> fluxScheme_;

        //- Cached advection changes
        tmp<volScalarField> alphaAdvection_;
        tmp<volScalarField> alphaRhoAdvection_;
        tmp<volVectorField> alphaRhoUAdvection_;
        tmp<volScalarField> alphaRhoEAdvection_;

        virtual tmp<volScalarField> ESource() const;


public:

    //- Runtime type information
    ClassName("fluid");

    // Constructor
    fluidPhaseModel
    (
        const phaseSystem& fluid,
        const word& phaseName,
        const label index,
        const word& thermoType = word::null
    );


    //- Destructor
    virtual ~fluidPhaseModel();


    // Member Functions

        // integrationSystem functions

            //- Decode primative variables
            virtual void decode();

            //- Encode conserved variables
            virtual void encode();

            //- Update fluxes
            virtual void update();

            //- Solve sub-step stepi
            virtual void solve();

            //- Solve viscous and source terms
            virtual void postExplicit();

            //- Solve viscous and source terms
            virtual void postImplicit();

            //- Store explicit deltas
            virtual void storeExplicit();

            //- Clear cached data
            virtual void clear();

            //- Correct delta terms
            virtual void correctDeltas();

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


        //- Return the density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Access the density
        virtual volScalarField& rho()
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

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return thermoPtr_->nu();
        }

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

        //- Return thermodynamic model
        virtual const ::Foam::fluidThermo& fluidThermo() const
        {
            return thermoPtr_();
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return thermoPtr_->nu(patchi);
        }

        //- Return the laminar viscosity for patch
        virtual scalar cellnu(const label celli) const
        {
            return thermoPtr_->cellnu(celli);
        }

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const
        {
            return thermoPtr_->mu();
        }

        //- Return the laminar dynamic viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return thermoPtr_->mu().boundaryField()[patchi];
        }

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const
        {
            return thermoPtr_->speedOfSound();
        }

        //- Return cell pressure
        virtual scalar cellpRhoT(const label celli) const
        {
            return thermoPtr_->cellpRhoT(celli);
        }

        //- Return cell speedOfSound
        virtual scalar cellSpeedOfSound
        (
            const scalar p,
            const label celli
        ) const
        {
            return thermoPtr_->cellSpeedOfSound(p, celli);
        }

        //- Correct thermodynamic quantities
        virtual void correctThermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
