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
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseCompressibleSystem
{
protected:
// Protected data

    //- Primitative variables

        //- Total mass
        volScalarField rho_;

        //- Velocity
        volVectorField U_;

        //- Internal energy
        volScalarField e_;

        //- Pressure
        volScalarField p_;


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


    // ODE variables

        //- Stored old field indexes
        labelList oldIs_;

        //- Stored delta indexes
        labelList deltaIs_;

        //- Old values for ode solver
        PtrList<volVectorField> rhoUOld_;
        PtrList<volScalarField> rhoEOld_;

        //- Stored delta fields
        PtrList<volVectorField> deltaRhoU_;
        PtrList<volScalarField> deltaRhoE_;

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

        //- Update fluxes
        virtual void update() = 0;

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


    // Member Access Functions

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

        //- Const access to internal energy
        const volScalarField& e() const
        {
            return e_;
        }

        //- Non-const access to internal energy
        volScalarField& e()
        {
            return e_;
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

        //- Return speed of sound
        virtual tmp<volScalarField> c() const = 0;

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
