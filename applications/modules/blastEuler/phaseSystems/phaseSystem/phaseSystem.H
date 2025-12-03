/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::phaseSystem

Description

SourceFiles
    phaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "PtrListDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "timeIntegrationSystem.H"
#include "timeIntegrator.H"
#include "phaseFluxScheme.H"
#include "HashPtrTable.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "masterSystemList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendedDragModel;
class blendedVirtualMassModel;
class blendedHeatTransferModel;
class blendedLiftModel;
class blendedWallLubricationModel;
class blendedTurbulentDispersionModel;
class massTransferModel;
class interfacialPressureModel;
class interfacialVelocityModel;
class pressureRelaxationModel;
class dragODE;
class pressureRelaxationSolver;

class blendingMethod;
template<class modelType>
class BlendedInterfacialModel;


/*---------------------------------------------------------------------------*\
                      Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary,
    public timeIntegrationSystem
{
public:

    // Public typedefs

        typedef PtrListDictionary<phaseModel> phaseModelList;

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef HashTable
        <
            autoPtr<blendedDragModel>,
            phasePairKey,
            phasePairKey::hash
        > dragModelTable;

        typedef HashTable
        <
            autoPtr<blendedVirtualMassModel>,
            phasePairKey,
            phasePairKey::hash
        > virtualMassModelTable;

        typedef HashTable
        <
            autoPtr<blendedLiftModel>,
            phasePairKey,
            phasePairKey::hash
        > liftModelTable;

        typedef HashTable
        <
            autoPtr<blendedWallLubricationModel>,
            phasePairKey,
            phasePairKey::hash
        > wallLubricationModelTable;

        typedef HashTable
        <
            autoPtr<blendedTurbulentDispersionModel>,
            phasePairKey,
            phasePairKey::hash
        > turbulentDispersionModelTable;

        typedef HashTable
        <
            autoPtr<blendedHeatTransferModel>,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;

        typedef HashTable
        <
            autoPtr<massTransferModel>,
            phasePairKey,
            phasePairKey::hash
        > massTransferModelTable;

        typedef HashTable
        <
            autoPtr<interfacialPressureModel>,
            phasePairKey,
            phasePairKey::hash
        > interfacialPressureModelTable;

        typedef HashTable
        <
            autoPtr<interfacialVelocityModel>,
            phasePairKey,
            phasePairKey::hash
        > interfacialVelocityModelTable;

        typedef HashTable
        <
            autoPtr<pressureRelaxationModel>,
            phasePairKey,
            phasePairKey::hash
        > pressureRelaxationModelTable;

        // Coefficient tables
        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > KdTable;

        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > VmTable;

        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > mDotTable;

        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > PiTable;

        typedef HashPtrTable
        <
            volVectorField,
            phasePairKey,
            phasePairKey::hash
        > UiTable;

        static const dimensionedScalar zeroMDot;

        enum PIPressure
        {
            VOLUME,
            TOTAL,
        };
        const static NamedEnum<PIPressure, 2> PIPressureNames_;

        enum PVRelaxation
        {
            NONE,
            MODEL,
            ODE,
            INSTANT
        };
        const static NamedEnum<PVRelaxation, 4> PVRelaxationNames_;


protected:

    // Protected typedefs

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        typedef
            HashTable<autoPtr<blendingMethod>, word, word::hash>
            blendingMethodTable;

        typedef
            HashTable
            <
                autoPtr<aspectRatioModel>,
                phasePairKey,
                phasePairKey::hash
            >
            aspectRatioModelTable;


    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture velocity
        volVectorField U_;

        //- Mixture flux
        surfaceScalarField phi_;

        //- Mixture mass flux
        surfaceScalarField rhoPhi_;

        //- Mixture pressure
        volScalarField p_;

        //- Interfacial pressure type
        PIPressure PIModel_;

        //- Total mixture pressure
        autoPtr<volScalarField> PIPtr_;

        //- Mixture temperature
        volScalarField T_;

        //- Mixture thermal diffusivity
        volScalarField kappa_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Fluid phase phase models
        UPtrList<phaseModel> fluidPhaseModels_;

        //- Slave phase phase models
        UPtrList<phaseModel> slavePhaseModels_;

        //- Pointer to kinetic theory system
        masterSystemList& master_;

        //- Phase pairs
        phasePairTable phasePairs_;

        //- Ode system and solver for solving the  drag system
        autoPtr<dragODE> dragODE_;

        //- Ode system and solver for solving the  drag system
        autoPtr<pressureRelaxationSolver> pressureSolver_;


        // Sub Models

            //- Velocity relaxation type
            PVRelaxation VRelaxation_;

            //- Pressure relaxation type
            PVRelaxation PRelaxation_;

            //- Blending methods
            blendingMethodTable blendingMethods_;

            //- Aspect ratio models
            aspectRatioModelTable aspectRatioModels_;

            //- Drag models
            dragModelTable dragModels_;

            //- Virtual mass models
            virtualMassModelTable virtualMassModels_;

            //- Lift models
            liftModelTable liftModels_;

            //- Wall lubrication models
            wallLubricationModelTable wallLubricationModels_;

            //- Turbulent dispersion models
            turbulentDispersionModelTable turbulentDispersionModels_;

            //- Heat transfer models
            heatTransferModelTable heatTransferModels_;

            //- Mass transfer models
            massTransferModelTable massTransferModels_;
	        List<List<bool>> hasMassTransfer_;

            //- Interfacial pressure models
            interfacialPressureModelTable interfacialPressureModels_;

            //- Interfacial velocity models
            interfacialVelocityModelTable interfacialVelocityModels_;

            //- Pressure relaxation models
            pressureRelaxationModelTable pressureRelaxationModels_;


        //- Stored coefficients

            //- Mass transfer rates
            mDotTable mDots_;


    // Private member functions

        //- Generate pairs
        void generatePairs
        (
            const dictTable& modelDicts
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void createSubModels
        (
            const dictTable& modelDicts,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool required
        );

        //- Generate pairs and blended sub-model tables
        template<class BlendedModel>
        void generateBlendedPairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<BlendedModel>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool required,
            const bool correctFixedFluxBCs = true
        );

        //- Generate pairs and two-sided sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                Pair<autoPtr<modelType>>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool required,
            const bool correctFixedFluxBCs = true
        );

        //- Solve finite pressure relaxation (not implemented yet)
//         void relaxPressure(const dimensionedScalar& deltaT);

        //- Solve finite velocity relaxation
        void relaxVelocity(const dimensionedScalar& deltaT);

        //- Solve heat transfer between phases
        void relaxTemperature(const dimensionedScalar& deltaT);

        //- Relax pressure between phases
        void relaxPressure(const dimensionedScalar& deltaT);


public:

    ClassName("phaseSystem");

    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh&);


    //- Destructor
    virtual ~phaseSystem();


    // Member Functions

        //- Calculate mixture properties
        void calcMixtureVariables();

        //- Calculate mixture fluxes
        void calcMixtureFluxes();

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Solve relaxation terms
        virtual void postExplicit();

        //- Solve viscous terms
        virtual void postImplicit();

        //- Store explicit changes
        virtual void storeExplicit();

        //- Post update
        virtual void postUpdate();

        //- Clear flux schemes
        virtual void clear();

        //- Print statistics
        virtual void printInfo() const;

        //- Return the mixture density
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the interfacial velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the interfacial pressure
        const volScalarField& p() const
        {
            return p_;
        }

        //- Return the interfacial pressure
        const volScalarField& PI() const
        {
            return PIPtr_.valid() ? PIPtr_() : p_;
        }

        //- Return the interfacial pressure
        scalar cellPI(const label celli) const;

        //- Return gravitational acceleration
        const uniformDimensionedVectorField& g() const
        {
            return g_;
        }

        //- Calculate and return the total volumetric flux
        tmp<surfaceScalarField> phi() const;

        //- Return the aspect-ratio for a pair
        tmp<volScalarField> E(const phasePairKey& key) const;

        //- Return the aspect-ratio for a pair
        scalar cellE(const phasePairKey& key, const label celli) const;

        //- Return is mas transfer is used for a phase
        bool hasMassTransfer(const phaseModel& phase) const;

        //- Return is mas transfer is used for a phase
        bool hasMassTransfer
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the mass transfer for a pair
        tmp<volScalarField> mDot
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the volume fraction source due to mass transfer for a pair
        tmp<volScalarField> mDotByRho
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the volume fraction source due to mass transfer for a pair
        tmp<volScalarField> mDotByRho
        (
            const volScalarField& mD,
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the momentum change due to mass transfer for a pair
        tmp<volVectorField> mDotU
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the momentum change due to mass transfer for a pair
        tmp<volVectorField> mDotU
        (
            const volScalarField& mD,
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the energy change due to mass transfer for a pair
        tmp<volScalarField> mDotE
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the energy change due to mass transfer for a pair
        tmp<volScalarField> mDotE
        (
            const volScalarField& mD,
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the pseudo thermal energy change due to
        //  mass transfer for a pair
        tmp<volScalarField> mDotPTE
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the pseudo thermal energy change due to
        //  mass transfer for a pair
        tmp<volScalarField> mDotPTE
        (
            const volScalarField& mD,
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Read base phaseProperties dictionary
        bool read();

        // Access

            // Sub-model lookup

            //- Return blending methods
            const blendingMethodTable& blendingMethods() const
            {
                return blendingMethods_;
            }

            //- Return a sub model between a phase pair
            template<class modelType>
            bool foundInterfacialModel(const phasePair& key) const;

            //- Return a sub model between a phase pair
            template<class modelType>
            const modelType& lookupInterfacialModel(const phasePair& key) const;

            //- Return a sub model between two phases
            template<class modelType>
            bool foundInterfacialModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous,
                const bool ordered = true
            ) const;

            //- Return a sub model between two phases
            template<class modelType>
            const modelType& lookupInterfacialModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous,
                const bool ordered = true
            ) const;

            //- Check availability of a blended sub model for a given phase pair
            template<class BlendedModelType>
            bool foundBlendedInterfacialModel(const phasePair& key) const;

            //- Return a blended sub model between a phase pair
            template<class BlendedModelType>
            const BlendedModelType&
            lookupBlendedInterfacialModel(const phasePair& key) const;

            //- Return the surface tension coefficient
            const dimensionedScalar& sigma() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return phase models
            inline const phaseModelList& phases() const;

            //- Return non-const access to phase models
            inline phaseModelList& phases();

            //- Return moving phase models
            inline const phaseModelList& movingPhases() const;

            //- Return the phase pairs
            inline const phasePairTable& phasePairs() const;

            //- Return the phase not given as an argument in a two-phase system
            //  An error is generated if the system is not two-phase
            inline const phaseModel& otherPhase(const phaseModel& phase) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
