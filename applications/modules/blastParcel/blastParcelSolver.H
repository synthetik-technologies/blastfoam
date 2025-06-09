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
    Foam::solvers::blastParcel

Description
    Module for blastParcelFoam solver

SourceFiles
    blastParcelSolver.C

\*---------------------------------------------------------------------------*/

#ifndef blastParcelSolver_H
#define blastParcelSolver_H

#include "explicitSolver.H"
#include "fvTimeIntegrator.H"
#include "coupledCompressibleSystem.H"
#include "parcelClouds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                        Class blastParcel Declaration
\*---------------------------------------------------------------------------*/

class blastParcel
:
    public explicitSolver
{
protected:
// Protected data


    //- Time integrator
    fvTimeIntegrator integrator_;

    //- Compressible system
    autoPtr<compressibleSystem> fluidPtr_;

    //- Clouds
    parcelClouds& clouds_;

    //- Volume fraction of the clouds
    volScalarField alphac_;

    // Protected member functions

        //- Read maximum Courant number
        virtual scalar CoNum() const;


public:


    //- Runtime type information
    TypeName("blastParcel");


    // Constructors

        //- Construct from region mesh
        blastParcel(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        blastParcel(const blastParcel&) = delete;


    //- Destructor
    virtual ~blastParcel();


    // Member Functions

        //- Solve
        virtual void solve();

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
