/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::regionSolvers::compressibleVof

Description
    Solver for a compressible phase using compressibleVofFoam libraries

SourceFiles
    compressibleVofRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleVofRegionSolver_H
#define compressibleVofRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "fluidRegionSolver.H"
#include "fvCFD.H"
#include "pimpleControl.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "noPhaseChange.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                           Class compressibleVofRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class compressibleVof
:
    public fluid
{
protected:
// Protected data

        //- Pimple
        pimpleControl pimple;

        //- Dynamic pressure
        volScalarField p_rgh;

        //- Velocity
        volVectorField U;

        //- Volumetric flux
        surfaceScalarField phi;

        //- Mixture thermo
        twoPhaseMixtureThermo mixture;

        //- Phase change model
        autoPtr<twoPhaseChangeModel> phaseChangePtr;
        twoPhaseChangeModel& phaseChange;

        //- Volume fraction2
        volScalarField& alpha1;
        volScalarField& alpha2;

        //- Phase densities
        const volScalarField& rho1;
        const volScalarField& rho2;

        //- Mixture density
        volScalarField& rho;

        //- Mixutre pressure
        volScalarField& p;

        //- Mixture temperature
        volScalarField& T;

        //- Phase compressibility
        const volScalarField& psi1;
        const volScalarField& psi2;

        // Pressure floor
        dimensionedScalar pMin;

        //- "readGravitationalAcceleration.H"
        //  "readhRef.H"
        //  "gh.H"
        uniformDimensionedVectorField g;
        uniformDimensionedScalarField hRef;
        dimensionedScalar ghRef;
        volScalarField gh;
        surfaceScalarField ghf;

        //- Pressure reference object
        Foam::pressureReference pressureReference;

        //- Mass flux
        surfaceScalarField rhoPhi;

        volScalarField dgdt;

        //- "createAlphaFluxes.H"
        bool alphaRestart;
        surfaceScalarField alphaPhi10;
        tmp<surfaceScalarField> talphaPhi1Corr0;
        tmp<surfaceScalarField> talphaPhi1;

        compressibleInterPhaseTransportModel turbulence;

        volScalarField K;

        mutable autoPtr<volScalarField> divU;
        mutable autoPtr<surfaceVectorField> Uf;
        mutable tmp<fvVectorMatrix> tUEqn;
        mutable tmp<volScalarField> rAU;
        mutable volScalarField::Internal contErr;

        bool correctPhi;
        scalar cumulativeContErr;

        IOMRFZoneList MRF;
        Foam::fvModels& fvModels;
        Foam::fvConstraints& fvConstraints;

        void solveU();
        void solvep();
        void solveT();
        void correctPhiField();


public:

    TypeName("compressibleVof");

    // Constructor
    compressibleVof(dynamicFvMesh& mesh, const regionSolverList& regions);

    //- Destructor
    virtual ~compressibleVof();


    // Member functions

        // Edit

            //- Move the mesh
            virtual bool moveMesh(const IterType iter);

            //- Solve the model
            virtual void solve();

            //- Return the Courant number
            virtual scalar CoNum() const;

            //- Return the maximum Courant number
            virtual scalar maxCo() const;

            //- Return the new time step
            virtual scalar newDeltaT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
