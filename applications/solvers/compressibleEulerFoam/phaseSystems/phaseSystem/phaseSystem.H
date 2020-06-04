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
#include "dragModel.H"
#include "integrationSystem.H"
#include "timeIntegrator.H"
#include "interfacialPressureModel.H"
#include "interfacialVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class virtualMassModel;
class heatTransferModel;
class liftModel;
class wallLubricationModel;
class turbulentDispersionModel;
class dragOde;
class pressureOde;
class kineticTheorySystem;

/*---------------------------------------------------------------------------*\
                      Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary,
    public integrationSystem
{
public:

    // Public typedefs

        typedef PtrListDictionary<phaseModel> phaseModelList;

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef HashTable
        <
            autoPtr<dragModel>,
            phasePairKey,
            phasePairKey::hash
        > dragModelTable;

        typedef HashTable
        <
            autoPtr<virtualMassModel>,
            phasePairKey,
            phasePairKey::hash
        > virtualMassModelTable;

        typedef HashTable
        <
            autoPtr<liftModel>,
            phasePairKey,
            phasePairKey::hash
        > liftModelTable;

        typedef HashTable
        <
            autoPtr<wallLubricationModel>,
            phasePairKey,
            phasePairKey::hash
        > wallLubricationModelTable;

        typedef HashTable
        <
            autoPtr<turbulentDispersionModel>,
            phasePairKey,
            phasePairKey::hash
        > turbulentDispersionModelTable;

        typedef HashTable
        <
            autoPtr<heatTransferModel>,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;


protected:

    // Protected typedefs

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

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

        //- Mixture pressure
        volScalarField p_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Interfacial pressure model
        autoPtr<interfacialPressureModel> interfacialPressure_;

        //- Interfacial velocity model
        autoPtr<interfacialVelocityModel> interfacialVelocity_;

        //- Pointer to kinetic theory system
        kineticTheorySystem* kineticTheoryPtr_;

        //- Is the kinetic theory polydisperse
        bool polydisperseKineticTheory_;

        //- Phase pairs
        phasePairTable phasePairs_;

        //- Ode system and solver for solving the  drag system
        autoPtr<dragOde> dragOde_;

        //- Switch to solve drag in ode solver
        Switch solveDragOde_;

        // Sub Models

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
            const bool
        );

        //- Solve finite pressure relaxation (not implemented yet)
//         void relaxPressure(const dimensionedScalar& deltaT);

        //- Solve finite velocity relaxation
        void relaxVelocity(const dimensionedScalar& deltaT);

        //- Solve heat transfer between phases
        void relaxTemperature(const dimensionedScalar& deltaT);

        //- Calculate the mixture density
        void calcRho();


public:

    ClassName("phaseSystem");

    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh&);


    //- Destructor
    virtual ~phaseSystem();


    // Member Functions

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

        //- Return the mixture density
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the mixture pressure
        const volScalarField& p() const
        {
            return p_;
        }

        //- Calculate and return the total volumetric flux
        tmp<surfaceScalarField> phi() const;

        //- Return the aspect-ratio for a pair
        tmp<volScalarField> E
        (
            const phasePairKey& key,
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Return the aspect-ratio for a pair
        scalar cellE
        (
            const label celli,
            const phasePairKey& key,
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Return the aspect-ratio for a pair
        tmp<volScalarField> E(const phasePairKey& key) const;

        //- Relax pressures and velocities
        void relax();

        //- Read base phaseProperties dictionary
        bool read();

        // Access

            // Sub-model lookup

            //- Return a sub model between a phase pair
            template<class modelType>
            bool foundSubModel(const phasePair& key) const;

            //- Return a sub model between a phase pair
            template<class modelType>
            const modelType& lookupSubModel(const phasePair& key) const;

            //- Return a sub model between two phases
            template<class modelType>
           bool foundSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;

            //- Return a sub model between two phases
            template<class modelType>
            const modelType& lookupSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;


            //- Return the surface tension coefficient
            const dimensionedScalar& sigma() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return phase models
            inline const phaseModelList& phases() const;

            //- Return non-const access to phase models
            inline phaseModelList& phases();

            //- Return the phase pairs
            inline const phasePairTable& phasePairs() const;
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
