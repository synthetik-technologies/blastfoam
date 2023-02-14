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
    Foam::regionSolvers::solid

Description
    Solver for a solid region solver

SourceFiles
    solidRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidRegionSolver_H
#define solidRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "solidModel.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                           Class regionSolver Declaration
\*---------------------------------------------------------------------------*/

class solid
:
    public regionSolver
{
protected:
// Protected data

        //- Solid model
        autoPtr<solidModel> solid_;


public:

    TypeName("solid");

    // Constructor
    solid(dynamicFvMesh& mesh, const regionSolverList& regions);

    //- Destructor
    virtual ~solid();


    // Member functions

        // Access

            //- Return if the region solver is solid
            inline bool isSolid() const
            {
                return true;
            }


        // Edit

            //- Initialise the mesh
            virtual void initialiseMesh(const IterType firstIter);

            //- Initialise fields
            virtual void initialiseFields();

            //- Initialise the solver
            virtual void initialise();

            //- Move the mesh
            virtual bool moveMesh(const IterType iter);

            //- Solve the model
            virtual void solve();

            //- Clear the solver
            virtual void clear(const bool full = false);

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
