/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::kineticTheorySystem

Description
    Class to handle the interactions of all granular phases. All models
    are stored within, and one model applies to all phases. The models are
    read from the kineticTheory subDict in phase properties.

SourceFiles
    kineticTheorySystem.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheorySystem_H
#define kineticTheorySystem_H

#include "phaseSystem.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "HashPtrTable.H"
#include "PtrList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class kineticTheoryModel;

namespace kineticTheoryModels
{
    class packingLimitModel;
    class radialModel;
    class viscosityModel;
    class frictionalStressModel;
    class granularPressureModel;
    class conductivityModel;
}

/*---------------------------------------------------------------------------*\
                 Class kineticTheorySystem Declaration
\*---------------------------------------------------------------------------*/

class kineticTheorySystem
:
    public regIOobject
{
    // Private data

        //- Const refrence to multiphaseSystem
        const phaseSystem& fluid_;

        //- Kinetic theory dictionary
        const dictionary& dict_;

        //- Name of total partical phase
        const word name_;

        //- Write the total granular phase
        Switch writeTotal_;

        // Total granualr fields

            //- Total granular phase-fraction
            volScalarField alphap_;

            //- Phase averaged velocity
            volVectorField Up_;

            //- Phase averaged granular temperature
            volScalarField Thetap_;


        //- PtrList of kinetic theory models
        PtrListDictionary<tmp<kineticTheoryModel>>
            kineticTheoryModels_;

        //- PtrList of granular temperatures
        PtrListDictionary<tmp<volScalarField>> Thetas_;

        //- List of granular phase indexes
        labelList phaseIndexes_;

        //- Table of granular pairs
        List<phasePairKey> pairs_;

        //- Run-time selected packing limit model
        autoPtr<kineticTheoryModels::packingLimitModel> packingLimitModel_;

        //- Run-time selected radial distribution model
        autoPtr<kineticTheoryModels::radialModel> radialModel_;

        //- Run-time selected viscosity model
        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        //- Run-time selected granular pressure model
        autoPtr<kineticTheoryModels::granularPressureModel>
            granularPressureModel_;

        //- Run-time selected thermal conductivity model
        autoPtr<kineticTheoryModels::conductivityModel>
            conductivityModel_;

        //- Run-time selected frictional stress model
        autoPtr<kineticTheoryModels::frictionalStressModel>
            frictionalStressModel_;


        // Kinetic Theory Model coefficients

            //- Table for coefficient of restitutions based on phase pairs
            phasePair::scalarTable eTable_;

            //- Table for coefficient of fricitions based on phase pairs
            phasePair::scalarTable CfTable_;

            //- Maximum packing phase-fraction
            volScalarField alphaMax_;

            //- Minimum phase alpha max
            scalar minAlphaMax_;

            //- Volume fraction that frictional stressed are turned on
            volScalarField alphaMinFriction_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;


        // Stored tables of coefficients

            //- Frictional pressure
            volScalarField Pfr_;

            //- Frictional pressure prime
            volScalarField PfrPrime_;

            //- Frictional viscosity
            volScalarField nuFric_;

            //- Minimum granular temperature
            dimensionedScalar minTheta_;


    // Private Member Functions

        //- Compute the packing limit for phase1 and phase2
        scalar calcAlphaMax(const label, const label, const scalarList&) const;


public:

    // Constructors

        //- Construct from components
        kineticTheorySystem(const phaseSystem&);


    //- Destructor
    ~kineticTheorySystem();

    //- Disallow default bitwise copy construct
    kineticTheorySystem(const kineticTheorySystem&) = delete;


    // Member Functions

        //- Return total particle phase name
        const word& name() const
        {
            return name_;
        }

        //- Return if granular system is polydiserse
        bool polydisperse() const;

        //- Re-read model coefficients if they have changed
        virtual bool read();

        virtual bool readIfModified();

        //- Return fluid
        const phaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Const refrence to list of granular indexes
        const labelList& phaseIndexes() const;

        //- Add phase to list of kinetic theory phases.
        void addPhase(kineticTheoryModel&);

        //- Return true if phase is granular
        bool found(const word&) const;

        //- Return residual alpha
        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        //- Return the minimum max alpha
        scalar minAlphaMax() const
        {
            return minAlphaMax_;
        }

        //- Return total granular phase volume fraction
        const volScalarField& alphap() const
        {
            return alphap_;
        }

        //- Return phase averaged velocity
        const volVectorField& Up() const
        {
            return Up_;
        }

        //- Return maximum volume fraction
        const volScalarField& alphaMax() const
        {
            return alphaMax_;
        }

        //- Return the minimum volume fraction frictional stresses are
        //  present
        const volScalarField& alphaMinFriction() const
        {
            return alphaMinFriction_;
        }

        //- Coefficients of restitutions between granular phases
        scalar es(const phasePairKey& pair) const;

        //- coefficients of friction between granular phases
        scalar Cf(const phasePairKey& pair) const;

        //- Return the list of granular temperatures
        const PtrList<tmp<volScalarField>>& Thetas() const
        {
            return Thetas_;
        }

        //- Return radial distribution function
        tmp<volScalarField> gs0
        (
            const phaseModel&,
            const phaseModel&,
            const bool calc = false
        ) const;

        //- Return radial distribution function
        tmp<volScalarField> gs0Prime
        (
            const phaseModel&,
            const phaseModel&,
            const bool calc = false
        ) const;

        //- Return the viscosity
        tmp<volScalarField> nu
        (
            const phaseModel&,
            const volScalarField&
        ) const;

        //- Return the sum of all granular pressure for a given phase
        tmp<volScalarField> Ps(const phaseModel&) const;

        //- Return the granular pressure for two phases
        tmp<volScalarField> Ps(const phaseModel&, const phaseModel&) const;

        //- Return the sum of the derivative of granular pressure
        //  for a given phase
        tmp<volScalarField> dPsdAlpha(const phaseModel&) const;

        //- Return the sum of the derivative of granular pressure
        //  for a given phase
        tmp<volScalarField> dPsdTheta(const phaseModel&) const;

        //- Return the thermal conductivity
        tmp<volScalarField> kappa
        (
            const phaseModel&,
            const volScalarField&
        ) const;

        //- Return the fricitional-pressure
        const volScalarField& frictionalPressure() const;

        //- Return the frictional-pressure'
        // (derivative of frictional-pressure w.r.t. phase-fraction)
        const volScalarField& frictionalPressurePrime
        (
            const volScalarField& alpha
        ) const;

        //- Return the frictional-viscosity'
        const volScalarField& nuFrictional() const;

        //- Solve the kinetic theory equations and correct the viscosity
        void correct();

        //- Correct volume fraction total
        void correctAlphap();

        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kineticTheorySystem&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //