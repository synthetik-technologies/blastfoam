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
    Foam::phaseCompressibleSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    phaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef psiuCompressibleSystem_H
#define psiuCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystem.H"
#include "fluxScheme.H"
#include "psiuMulticomponentThermo.H"
#include "RASThermophysicalTransportModel.H"
#include "unityLewisEddyDiffusivity.H"
#include "uniformDimensionedFields.H"
#include "laminarFlameSpeed.H"
#include "ignition.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class psiuCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class psiuCompressibleSystem
:
    public compressibleSystem
{
protected:
// Protected data

    //- Thermo
    autoPtr<psiuMulticomponentThermo> thermo_;

    //- Primitative variables

        //- Total mass
        volScalarField rho_;

        //- Pressure
        volScalarField& p_;

        //- Temperature
        volScalarField& T_;

        //- Internal energy
        volScalarField& e_;

        //- Refernce to unburnt internal energy
        volScalarField& eu_;

    //- Conserved variables

        //- Total unbunrt energy
        volScalarField rhoEu_;


    // - Flame models

        //- Combustion properties
        IOdictionary combustionProperties_;

        //- Burnt fraction
        volScalarField& b_;

        //- Velocity ratio
        volScalarField Xi_;

        //- Strained laminar flame speed
        volScalarField Su_;

        //- Turbulent flame speed
        volScalarField St_;

        //- Su bounds
        dimensionedScalar SuMin_;
        dimensionedScalar SuMax_;

        // - Laminar flame speed model
        autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed_;

        // Model names
        word SuModel_;
        word XiModel_;

        dimensionedScalar sigmaExt_;
        dimensionedScalar XiCoef_;
        dimensionedScalar XiShapeCoef_;
        dimensionedScalar uPrimeCoef_;

        //- Ignition model
        ignition ign_;

        // Thermophysical transport
        autoPtr
        <
            turbulenceThermophysicalTransportModels::unityLewisEddyDiffusivity
            <
                RASThermophysicalTransportModel
                <
                    ThermophysicalTransportModel
                    <
                        compressibleMomentumTransportModel,
                        fluidThermo
                    >
                >
            >
        > thermophysicalTransport_;


        //- Calculate and return the turbulent flame-speed kernel correction
        dimensionedScalar calcStCorr
        (
            const volScalarField& c,
            const surfaceScalarField& nf,
            const dimensionedScalar& dMgb
        ) const;


public:

    TypeName("psiuCompressibleSystem");

    // Constructor
    psiuCompressibleSystem
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~psiuCompressibleSystem();


    // Member functions

        //- Non-const reference to thermo
        fluidThermo& thermo()
        {
            return thermo_();
        }

        //- Const reference to thermo
        const fluidThermo& thermo() const
        {
            return thermo_();
        }

        //- Non-const reference to reaction thermo
        psiuMulticomponentThermo& rThermo()
        {
            return thermo_();
        }

        //- Const reference to reaction thermo
        const psiuMulticomponentThermo& rThermo() const
        {
            return thermo_();
        }

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Set old lists and fluxes (initialization of fields)
        virtual void postUpdate();


    // Member Access Functions

        //- Const access to total mass
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Access to total mass
        volScalarField& rho()
        {
            return rho_;
        }

        //- Const access to total mass
        const volScalarField& rhoEff() const
        {
            return rho_;
        }

        //- Access to total mass
        volScalarField& rhoEff()
        {
            return rho_;
        }

        //- Const access to internal energy
        virtual const volScalarField& he() const
        {
            return e_;
        }

        //- Non-const access to internal energy
        virtual volScalarField& he()
        {
            return e_;
        }

        //- Return thermodynamic pressure
        const volScalarField& p() const
        {
            return p_;
        }

        //- Return temperature
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
