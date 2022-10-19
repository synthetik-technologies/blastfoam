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
    Foam::compressibleSystem

Description
    Base class for a region solver

SourceFiles
    regionSolver.C
    regionSolverNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionSolver_H
#define regionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dynamicFvMesh.H"
#include "globalPolyBoundaryMesh.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regionSolver Declaration
\*---------------------------------------------------------------------------*/

class regionSolver
{
protected:
// Protected data

        // Time
        const Time& runTime_;

        //- Dynamic mesh
        dynamicFvMesh& dynMesh_;

        //- fvMesh
        fvMesh& mesh_;

        //- Global polyBoundaryMesh
        globalPolyBoundaryMesh& globalBoundary_;

public:

    TypeName("regionSolver");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            regionSolver,
            dictionary,
            (
                dynamicFvMesh& mesh
            ),
            (mesh)
        );

    // Constructor
    regionSolver(dynamicFvMesh& mesh);


    // Selector
        static autoPtr<regionSolver> New
        (
            const word& type,
            dynamicFvMesh& mesh
        );


    //- Destructor
    virtual ~regionSolver();


    // Member functions

        // Access Functions

            //- Return the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Access the mesh
            inline fvMesh& mesh()
            {
                return mesh_;
            }


        // Edit

            //- Update global patches
            virtual void update();

            //- Apply topological changes
            virtual bool changeMesh();

            //- Initialise the mesh
            virtual void initialiseMesh() = 0;

            //- Initialise the solver
            virtual void initialise() = 0;

            //- Move the mesh
            virtual bool moveMesh(const bool finalIter);

            //- Solve the model
            virtual void solve() = 0;

            //- Return the Courant number
            virtual scalar CoNum() const = 0;

            //- Return the maximum Courant number
            virtual scalar maxCo() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
