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
    Foam::compressibleSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    compressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSystem_H
#define compressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxScheme.H"
#include "timeIntegrationSystem.H"
#include "runTimeSelectionTables.H"
#include "fluidThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "uniformDimensionedFields.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "hashedWordList.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class compressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class compressibleSystem
:
    public timeIntegrationSystem
{
protected:
// Protected data

        //- Velocity
        volVectorField U_;

    //- Conserved variables

        //- Momentum
        volVectorField rhoU_;

        //- Energy
        volScalarField rhoE_;


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

        //- Thermophysical transport model
        autoPtr<fluidThermophysicalTransportModel>
            thermophysicalTransport_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;

        //- Solution directions
        vector solutionDs_;

        //- Allocate models
        virtual void setModels();

public:

    TypeName("compressibleSystem");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            compressibleSystem,
            singlePhase,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            compressibleSystem,
            twoPhase,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            compressibleSystem,
            multiphase,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );

    // Constructor
    compressibleSystem(const fvMesh& mesh);


    // Selector
        static autoPtr<compressibleSystem> New(const fvMesh& mesh);

        template<class SelectionTable>
        static autoPtr<compressibleSystem> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& defaultType,
            SelectionTable* table
        );

    //- Destructor
    virtual ~compressibleSystem();


    // Member functions

        //- Return the flux scheme
        inline const fluxScheme& flux() const;

    // Member Access Functions

        //- Return if the system is inviscid
        inline bool inviscid() const;

        //- Return turbulence class
        inline const compressible::momentumTransportModel& turbulence() const;

        //- Return turbulence class
        inline compressible::momentumTransportModel& turbulence();

        //- Const access to fluidThermophysicalTransportModel model
        inline const fluidThermophysicalTransportModel&
        thermophysicalTransport() const;

        //- Non-const access to fluidThermophysicalTransportModel model
        inline fluidThermophysicalTransportModel& thermophysicalTransport();

        //- Const access to thermo model
        virtual const fluidThermo& thermo() const = 0;

        //- Non-const access to thermo model
        virtual fluidThermo& thermo() = 0;

        //- Const access to total mass
        virtual const volScalarField& rho() const = 0;

        //- Const access to velocity
        inline const volVectorField& U() const;

        //- Const access to internal energy
        virtual const volScalarField& he() const = 0;

        //- Access to internal energy
        virtual volScalarField& he() = 0;

        //- Const access to pressure
        virtual const volScalarField& p() const = 0;

        //- Return calculated temperature field
        virtual const volScalarField& T() const = 0;

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const = 0;

        //- Return volumetric flux
        inline const surfaceScalarField& phi() const;

        //- Return total mass flux
        inline const surfaceScalarField& rhoPhi() const;


    // Edit

        //- Decode primative variables
        virtual void decode() = 0;

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Clear the flux scheme
        virtual void clear();

        //- Return the maximum Courant number
        virtual scalar CoNum() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystemI.H"

#ifdef NoRepository
#include "compressibleSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
