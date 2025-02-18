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

    Reference:
    \verbatim
        Huilin, L., Gidaspow, D. (2003).
        Hydrodynamics of binary fluidization in a riser: CFD simulation using two
        granular temperatures.
        Chemical Engineering Science 58, 3777–3792.
    \endverbatim

SourceFiles
    kineticTheorySystem.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheorySystem_H
#define kineticTheorySystem_H

#include "masterSystem.H"
#include "phaseSystem.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "HashPtrTable.H"
#include "UPtrList.H"


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
    public masterSystem
{
    // Private data

        //- Phase averaged granular temperature
        autoPtr<volScalarField> ThetapPtr_;

        //- PtrList of kinetic theory models
        UPtrList<kineticTheoryModel> kineticTheoryModels_;

        //- PtrList of granular temperatures
        UPtrList<volScalarField> Thetas_;

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


        // Stored tables of coefficients

            //- Minimum granular temperature
            dimensionedScalar minTheta_;

            //- Are viscous effects included
            bool includeViscosity_;


    // Private Member Functions

        //- Compute the packing limit for phase1 and phase2
        scalar calcAlphaMax(const label, const label, const scalarList&) const;


public:

    TypeName("kineticTheorySystem");

    // Constructors

        //- Construct from components
        kineticTheorySystem(const phaseSystem&);

    //- Destructor
    ~kineticTheorySystem();

    //- Disallow default bitwise copy construct
    kineticTheorySystem(const kineticTheorySystem&) = delete;


    // Member Functions

        //- Return true if viscosity is included
        bool includeViscosity() const
        {
            return includeViscosity_;
        }

        //- Add phase to list of kinetic theory phases.
        virtual void addPhase(phaseModel&);

        //- Return the minimum max alpha
        scalar minAlphaMax() const
        {
            return minAlphaMax_;
        }

        //- Return maximum volume fraction
        virtual tmp<volScalarField> alphaMax() const
        {
            return alphaMax_;
        }

        //- Return maximum volume fraction
        virtual scalar alphaMax(const label celli) const
        {
            return alphaMax_[celli];
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
        const UPtrList<volScalarField>& Thetas() const
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

        //- Return radial distribution function for a cell
        scalar cellgs0
        (
            const label,
            const phaseModel&,
            const phaseModel&
        ) const;

        //- Return radial distribution function
        tmp<volScalarField> gs0Prime
        (
            const phaseModel&,
            const phaseModel&,
            const bool calc = false
        ) const;

        //- Return radial distribution function for a cell
        scalar cellgs0Prime
        (
            const label,
            const phaseModel&,
            const phaseModel&
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

        //- Return the bulk viscosity for a given phase
        tmp<volScalarField> lambda(const phaseModel& phase) const;

        //- Return the fricitional-pressure
        tmp<volScalarField> frictionalPressure(const phaseModel&) const;

        //- Return the frictional-pressure'
        // (derivative of frictional-pressure w.r.t. phase-fraction)
        tmp<volScalarField> frictionalPressurePrime
        (
            const phaseModel&
        ) const;

        //- Return the frictional-viscosity'
        tmp<volScalarField> muFrictional
        (
            const phaseModel&,
            const volScalarField&
        ) const;

        //- Return the PTE production
        tmp<volScalarField> productionSource
        (
            const kineticTheoryModel&,
            const phaseModel&
        ) const;

        //- Return the PTE dissipation rate
        tmp<volScalarField> dissipationSource
        (
            const phaseModel&,
            const phaseModel&,
            const dimensionedScalar& deltaT
        ) const;

        //- Initialize the model after everything is constructed
        virtual void initialize();

        //- Update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Remove stored fields
        virtual void clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kineticTheorySystem&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
