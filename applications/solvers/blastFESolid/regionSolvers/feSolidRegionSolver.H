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

#ifndef feSolidRegionSolver_H
#define feSolidRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"

#include "feMesh1.H"
#include "pointFields.H"
#include "materialModel.H"
#include "displacementConstraint.H"
#include "boundaryTraction.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                           Class regionSolver Declaration
\*---------------------------------------------------------------------------*/

class feSolid
:
    public regionSolver
{
protected:
// Protected data

    //- Finite element mesh
    feMesh1 femesh_;

    //- Solid properties
    IOdictionary solidProperties_;

    //- Mechanical properties
    IOdictionary mechanicalProperties_;

    materialModel::GeoType geoType_;

    //- Displacement field
    pointVectorField D_;

    //- Velocity field
    pointVectorField U_;

    //- Characteristic length of elements
    scalarField L_;

    //- Force
    pointVectorField force_;

    //- Mass matrix (row summed)
    pointScalarField::Internal M_;

    //- Cauchy stress tensor
    pointSymmTensorField sigma_;

    //- pressure field
    pointScalarField pressure_;

    //- Maximum wavespeed
    scalar maxVp_;

    //- Materials
    PtrList<materialModel> materials_;

    //- Displacement constraints
    PtrList<displacementConstraint> constraints_;

    //- Boundary loads
    PtrList<boundaryTraction> loads_;



public:

    TypeName("feSolid");

    // Constructor
    feSolid(fvMesh& mesh, const regionSolverList& regions);

    //- Destructor
    virtual ~feSolid();


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

            //- Return the new deltaT
            virtual scalar newDeltaT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
