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
    Foam::solvers::blastSolid

Description
    Module for blastSolidFoam solver

SourceFiles
    blastSolidSolver.C

\*---------------------------------------------------------------------------*/

#ifndef blastSolidSolver_H
#define blastSolidSolver_H

#include "explicitSolver.H"
#include "fvTimeIntegrator.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                        Class blastSolid Declaration
\*---------------------------------------------------------------------------*/

class blastSolid
:
    public explicitSolver
{
protected:
// Protected data

    //- Compressible system
    autoPtr<solidModel> solidPtr_;

    // Protected member functions

        //- Read maximum Courant number
        virtual scalar CoNum() const;

        //- Read maximum Diffusion number
        virtual scalar DiNum() const;


public:


    //- Runtime type information
    TypeName("blastSolid");


    // Constructors

        //- Construct from region mesh
        blastSolid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        blastSolid(const blastSolid&) = delete;


    //- Destructor
    virtual ~blastSolid();


    // Member Functions

        //- Solve
        virtual void solveExplicit();
        virtual void solveImplicit();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
