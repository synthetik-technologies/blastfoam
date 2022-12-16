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
    Foam::regionSolvers::thermal

Description
    Base class for a thermal region solver

SourceFiles
    thermalRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef thermalRegionSolver_H
#define thermalRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "solidBlastThermo.H"
#include "fvConstraints.H"
#include "fvModels.H"


namespace Foam
{
namespace regionSolvers
{
/*---------------------------------------------------------------------------*\
                           Class thermalRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class thermal
:
    public regionSolver
{
protected:
// Protected data

    //- Solid thermo model
    autoPtr<solidBlastThermo> thermo_;

    //- Fv models
    fvModels& fvModels_;

    //- Fv constraints
    fvConstraints& fvConstraints_;


public:

    TypeName("thermal");

    // Constructor
    thermal(dynamicFvMesh& mesh, const regionSolverList& regions);


    //- Destructor
    virtual ~thermal();


    // Member functions

        // Edit

            //- Return the region is solid
            virtual bool isSolid() const
            {
                return true;
            }

            //- Initialise the mesh
            virtual void initialiseMesh(const IterType firstIter);

            //- Initialise the solver
            virtual void initialise();

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
