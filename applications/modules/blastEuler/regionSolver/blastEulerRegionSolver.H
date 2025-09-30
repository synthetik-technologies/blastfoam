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
    Foam::regionSolvers::blastEuler

Description
    Solver for a reacting fluid phase using standard OpenFOAM thermodynamics

SourceFiles
    blastEulerRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef blastEulerRegionSolver_H
#define blastEulerRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidRegionSolver.H"
#include "fvTimeIntegrator.H"
#include "phaseSystem.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                    Class reactingRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class blastEuler
:
    public fluid
{
    //- Time integrator
    fvTimeIntegrator integrator_;

    //- Compressible system
    phaseSystem fluid_;


public:

    TypeName("blastEuler");

    // Constructor
    blastEuler
    (
        fvMesh& mesh,
        const regionSolverList& regions
    );

    //- Destructor
    virtual ~blastEuler();


    // Member functions

        // Edit

            //- Apply topological changes
            virtual bool changeMesh();

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
