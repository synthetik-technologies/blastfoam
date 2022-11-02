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
    Foam::regionSolvers::incompressibleFluid

Description
    Solver for a compressible phase using incompressibleFluidFoam libraries

SourceFiles
    incompressibleFluidRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleFluidRegionSolver_H
#define incompressibleFluidRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "fluidRegionSolver.H"
#include "fvCFD.H"
#include "pimpleControl.H"
#include "singlePhaseTransportModel.H"
#include "uniformDimensionedFields.H"
#include "kinematicMomentumTransportModel.H"
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
                           Class incompressibleFluidRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class incompressibleFluid
:
    public fluid
{
protected:
// Protected data

        pimpleControl pimple_;

        //- Gravity
        uniformDimensionedVectorField g_;

        volVectorField U;
        surfaceScalarField phi;
        volScalarField p;
        pressureReference pressureReference_;
        singlePhaseTransportModel laminarTransport_;
        autoPtr<incompressible::momentumTransportModel> turbulence_;

        dimensionedScalar initialMass;

        autoPtr<volScalarField> divU;
        autoPtr<surfaceVectorField> Uf;
        mutable tmp<fvVectorMatrix> tUEqn;

        scalar cumulativeContErr;

        IOMRFZoneList MRF;
        Foam::fvModels& fvModels;
        Foam::fvConstraints& fvConstraints;

        void solveU();
        void solvep();
        void correctPhi();


public:

    TypeName("incompressibleFluid");

    // Constructor
    incompressibleFluid(dynamicFvMesh& mesh);

    //- Destructor
    virtual ~incompressibleFluid();


    // Member functions

        // Edit

            //- Move the mesh
            virtual bool move(const bool finalIter);

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
