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
    Foam::regionSolvers::fluid

Description
    Base class for a fluid region solver

SourceFiles
    fluidRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fluidRegionSolver_H
#define fluidRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "displacementRelaxation.H"
#include "volFields.H"
#include "pointFields.H"

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                           Class fluidRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class fluid
:
    public regionSolver
{
protected:
// Protected data

        //- Viscous stresses
        volSymmTensorField devRhoReff_;

        //- old point field
        autoPtr<pointField> pointsOldPtr_;

        //- Cell velocity
        autoPtr<volVectorField> cellDPtr_;

        //- Point velocity
        autoPtr<pointVectorField> pointDPtr_;

        //- Motion diffusivity
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Relaxation
        autoPtr<displacementRelaxation> relaxation_;

        //- Velocity field names
        wordList velocityFields_;

        //- Tolerances
        scalar tolerance_;
        scalar relTol_;


public:

    TypeName("fluid");

    // Constructor
    fluid(dynamicFvMesh& mesh, const regionSolverList& regions);


    //- Destructor
    virtual ~fluid();


    // Member functions

        // Edit

            //- Return if the solver converged
            virtual Convergence converged() const
            {
                if (tolerance_ < 0)
                {
                    return UNKNOWN_CONVERGENCE;
                }
                else if
                (
                    relaxation_->error() < tolerance_
                 || relaxation_->relError() < relTol_
                )
                {
                    return CONVERGED;
                }
                return NOT_CONVERGED;
            }

            //- Initialise the mesh
            virtual void initialiseMesh(const IterType iter);

            //- Initialise the solver
            virtual void initialise();

            //- Apply topological changes
            virtual bool changeMesh();

            //- Move the mesh
            virtual bool moveMesh(const IterType iter);

            //- Clear the model
            virtual void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
