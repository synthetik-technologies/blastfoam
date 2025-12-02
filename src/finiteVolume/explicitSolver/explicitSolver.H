/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::explicitSolver

Description
    Base class for explicit solvers

SourceFiles
    explicitSolver.C

\*---------------------------------------------------------------------------*/

#ifndef explicitSolver_H
#define explicitSolver_H

#include "solver.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                        Class explicitSolver Declaration
\*---------------------------------------------------------------------------*/

class explicitSolver
:
    public solver
{
protected:
// Protected member data

    // Control parameters

        //- Maximum allowed Courant number
        scalar maxCo_;

        //- Maximum allowed Diffusion number
        scalar maxDi_;

        //- Maximum time-step
        scalar maxDeltaT_;

        //- Gravitational acceleration
        uniformDimensionedVectorField g_;

protected:

    // Protected member functions

        //- Read controls
        virtual bool read();

        //- Return the current Courant number
        virtual scalar CoNum() const = 0;

        //- Return the current Diffusion number
        virtual scalar DiNum() const = 0;

public:


    //- Runtime type information
    TypeName("explicitSolver");


    // Constructors

        //- Construct from region mesh
        explicitSolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        explicitSolver(const explicitSolver&) = delete;


    //- Destructor
    virtual ~explicitSolver();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Solve explicit
        virtual void solveExplicit() = 0;

        //- Solve Implicit
        virtual void solveImplicit() = 0;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the thermophysical property equations,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
