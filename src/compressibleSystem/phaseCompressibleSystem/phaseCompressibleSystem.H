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

#ifndef phaseCompressibleSystem_H
#define phaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxScheme.H"
#include "integrationSystem.H"
#include "runTimeSelectionTables.H"
#include "blastCompressibleTurbulenceModelFwd.H"
#include "radiationModel.H"
#include "uniformDimensionedFields.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseCompressibleSystem
:
    public integrationSystem
{
protected:
// Protected data

    //- Primitative variables

        //- Total mass
        volScalarField rho_;

        //- Velocity
        volVectorField U_;

        //- Pressure
        volScalarField p_;

        //- Temperature
        volScalarField T_;

        //- Internal energy
        volScalarField e_;


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
        autoPtr<blastCompressibleTurbulenceModel> turbulence_;

        //- Radiation model
        autoPtr<radiationModel> radiation_;


    // ODE variables

        //- Stored old field indexes
        labelList oldIs_;
        label nOld_;

        //- Stored delta indexes
        labelList deltaIs_;
        label nDelta_;

        //- Old values for ode solver
        PtrList<volVectorField> rhoUOld_;
        PtrList<volScalarField> rhoEOld_;

        //- Stored delta fields
        PtrList<volVectorField> deltaRhoU_;
        PtrList<volScalarField> deltaRhoE_;

        //- External energy source
        tmp<volScalarField::Internal> extESource_;

        //- Drag matrix
        tmp<fvVectorMatrix> dragSource_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;

        //- Minimum temperature
        dimensionedScalar TLow_;


        //- Allocate turbulence model if mu >= 0
        virtual void setModels(const dictionary& dict);

        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho() = 0;

public:

    TypeName("phaseCompressibleSystem");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseCompressibleSystem,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructor
    phaseCompressibleSystem(const fvMesh& mesh, const dictionary& dict);


    // Selector
        static autoPtr<phaseCompressibleSystem> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~phaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode() = 0;

        //- Encode conserved variables
        virtual void encode() = 0;

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const = 0;

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

        //- Add external energy source
        void addESource(const volScalarField::Internal& extEsrc);

        //- Add to U coefficient
        void addUCoeff(const volScalarField::Internal& coeff);

        //- Add to U source
        void addUSource(const volVectorField::Internal& src);


    // Member Access Functions

        //- Return turbulence class
        const blastCompressibleTurbulenceModel& turbulence() const;

        //- Return turbulence class
        blastCompressibleTurbulenceModel& turbulence();

        //- Const access to thermo model
        virtual const fluidThermoModel& thermo() const = 0;

        //- Non-const access to thermo model
        virtual fluidThermoModel& thermo() = 0;

        //- Const access to total mass
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Non-const access to total mass
        volScalarField& rho()
        {
            return rho_;
        }

        //- Const access to velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Non-const access to velocity
        volVectorField& U()
        {
            return U_;
        }

        //- Const access to pressure
        const volScalarField& p() const
        {
            return p_;
        }

        //- Non-const access to pressure
        volScalarField& p()
        {
            return p_;
        }

        //- Return calculated temperature field
        virtual tmp<volScalarField> T() const
        {
            return T_;
        }

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
        virtual tmp<volScalarField> speedOfSound() const = 0;

        //- Const access to momentum
        const volVectorField& rhoU() const
        {
            return rhoU_;
        }

        //- Non-const access to momentum
        volVectorField& rhoU()
        {
            return rhoU_;
        }

        //- Const access to total energy
        const volScalarField& rhoE() const
        {
            return rhoE_;
        }

        //- Non-const access to total energy
        volScalarField& rhoE()
        {
            return rhoE_;
        }

        //- Return volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return total mass flux
        const surfaceScalarField& rhoPhi() const
        {
            return rhoPhi_;
        }

        //- Return total momentum flux
        const surfaceVectorField& rhoUPhi() const
        {
            return rhoUPhi_;
        }

        //- Return total energy flux
        const surfaceScalarField& rhoEPhi() const
        {
            return rhoEPhi_;
        }

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
