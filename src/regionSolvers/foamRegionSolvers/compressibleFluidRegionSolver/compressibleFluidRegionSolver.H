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
    Foam::regionSolvers::compressibleFluid

Description
    Solver for a compressible phase using compressibleFluidFoam libraries

SourceFiles
    compressibleFluidRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleFluidRegionSolver_H
#define compressibleFluidRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "fluidRegionSolver.H"
#include "fvCFD.H"
#include "pimpleControl.H"
#include "fluidThermo.H"
#include "uniformDimensionedFields.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                           Class compressibleFluidRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class compressibleFluid
:
    public fluid
{
protected:
// Protected data

        pimpleControl pimple_;

        //- Gravity
        uniformDimensionedVectorField g_;

        autoPtr<fluidThermo> thermo_;
        volScalarField rho;
        volVectorField U;
        surfaceScalarField phi;
        volScalarField& p;
        const volScalarField& psi;
        pressureReference pressureReference_;
        autoPtr<compressible::momentumTransportModel> turbulence_;
        autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport_;
        volScalarField dpdt;
        volScalarField K;

        dimensionedScalar initialMass_;

        autoPtr<volScalarField> divRhoU_;
        autoPtr<surfaceVectorField> rhoUf;
        mutable tmp<fvVectorMatrix> tUEqn;

        scalar cumulativeContErr;

        IOMRFZoneList MRF;
        Foam::fvModels& fvModels;
        Foam::fvConstraints& fvConstraints;

        void solveU();
        void solvep();
        void solveE();
        void correctPhi();


public:

    TypeName("compressibleFluid");

    // Constructor
    compressibleFluid(dynamicFvMesh& mesh, const regionSolverList& regions);

    //- Destructor
    virtual ~compressibleFluid();


    // Member functions

        // Edit

            //- Solve the model
            virtual void solve();

            //- Return the Courant number
            virtual scalar CoNum() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
