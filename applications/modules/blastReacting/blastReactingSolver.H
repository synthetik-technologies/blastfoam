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
    Foam::solvers::blastReacting

Description
    Module for blastReactingFoam solver

SourceFiles
    blastReactingSolver.C

\*---------------------------------------------------------------------------*/

#ifndef blastReactingSolver_H
#define blastReactingSolver_H

#include "explicitSolver.H"
#include "fvTimeIntegrator.H"
#include "reactingCompressibleSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                        Class blastReacting Declaration
\*---------------------------------------------------------------------------*/

class blastReacting
:
    public explicitSolver
{
protected:
// Protected data


    //- Time integrator
    fvTimeIntegrator integrator_;

    //- Compressible system
    reactingCompressibleSystem fluid_;

    // Protected member functions

        //- Read maximum Courant number
        virtual scalar CoNum() const;


public:


    //- Runtime type information
    TypeName("blastReacting");


    // Constructors

        //- Construct from region mesh
        blastReacting(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        blastReacting(const blastReacting&) = delete;


    //- Destructor
    virtual ~blastReacting();


    // Member Functions

        //- Solve
        virtual void solve();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
