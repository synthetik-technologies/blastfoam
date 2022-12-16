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
    Foam::regionSolvers::vof

Description
    Solver for a compressible phase using vofFoam libraries

SourceFiles
    vofRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef vofRegionSolver_H
#define vofRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "fluidRegionSolver.H"
#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "noPhaseChange.H"
#include "kinematicMomentumTransportModel.H"
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
                           Class vofRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class vof
:
    public fluid
{
protected:
// Protected data

        pimpleControl pimple;

        volScalarField p_rgh;
        volVectorField U;
        surfaceScalarField phi;

        immiscibleIncompressibleTwoPhaseMixture mixture;
        // autoPtr<twoPhaseChangeModel> phaseChangePtr;
        // twoPhaseChangeModel& phaseChange;
        volScalarField& alpha1;
        volScalarField& alpha2;

        const dimensionedScalar& rho1;
        const dimensionedScalar& rho2;
        volScalarField rho;

        bool alphaRestart;

        surfaceScalarField alphaPhi10;
        tmp<surfaceScalarField> talphaPhi1Corr0;
        surfaceScalarField rhoPhi;

        autoPtr<incompressible::momentumTransportModel> turbulence;

        //- Gravity
        uniformDimensionedVectorField g;
        uniformDimensionedScalarField hRef;
        dimensionedScalar ghRef;
        volScalarField gh;
        surfaceScalarField ghf;

        volScalarField p;

        Foam::pressureReference pressureReference;

        mutable autoPtr<volScalarField> divU;
        mutable autoPtr<surfaceVectorField> Uf;
        mutable tmp<fvVectorMatrix> tUEqn;
        mutable tmp<volScalarField> rAU;

        bool correctPhi;
        scalar cumulativeContErr;

        IOMRFZoneList MRF;
        Foam::fvModels& fvModels;
        Foam::fvConstraints& fvConstraints;

        void solveU();
        void solvep();
        void correctPhiField();


public:

    TypeName("vof");

    // Constructor
    vof(dynamicFvMesh& mesh, const regionSolverList& regions);

    //- Destructor
    virtual ~vof();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
