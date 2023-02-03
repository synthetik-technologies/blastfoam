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
    Foam::regionSolverList

Description
    List of region solvers

SourceFiles
    regionSolverList.C

\*---------------------------------------------------------------------------*/

#ifndef regionSolverList_H
#define regionSolverList_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListDictionary.H"
#include "Time.H"
#include "dynamicFvMesh.H"


namespace Foam
{
class regionSolver;

/*---------------------------------------------------------------------------*\
                           Class regionSolverList Declaration
\*---------------------------------------------------------------------------*/

// public enumerations
enum IterType
{
    FIRST_ITER,
    MID_ITER,
    FINAL_ITER
};

enum Convergence
{
    NOT_CONVERGED,
    CONVERGED,
    UNKNOWN_CONVERGENCE
};

class regionSolverList
:
    public PtrListDictionary<regionSolver>
{
protected:
// Protected data

        // Time
        const Time& runTime_;

        const dictionary& regionProperties_;

        //- Have the meshes changed
        boolList changed_;

        //- Predict solids
        bool predictSolids_;

        //- Is mapping fixed
        //  topological change require updating mapping
        bool fixedMapping_;

        //- Corrector index
        label iterNo_;

        //- Solve order
        Map<word> order_;

        //Has any region changed
        bool anyChanged() const
        {
            forAll(changed_, i)
            {
                if (changed_[i])
                {
                    return true;
                }
            }
            return false;
        }


public:

    TypeName("regionSolverList");

    // Constructor
    regionSolverList
    (
        const dictionary& regionProperties,
        const Time& runTime,
        PtrList<dynamicFvMesh>& regionMeshes,
        const List<Pair<word>>& regionTypes
    );


    //- Destructor
    virtual ~regionSolverList();


    // Member functions

        // Access

            //- Return solution controls dictionary
            const dictionary& solutionControls() const
            {
                return regionProperties_.subDict("solutionControls");
            }

            //- Return if all regions are converged
            bool converged() const;

            //- Return convergence type
            Convergence convergence() const;

            //- Return the Courant number
            virtual scalar CoNum() const;

            //- Return the maximum Courant number
            virtual scalar maxCo() const;

            //- Return the new deltaT
            virtual scalar newDeltaT() const;

            //- Return the current iteration index
            label iterNo() const
            {
                return iterNo_;
            }


        // Edit

            //- Initial update to fields
            virtual void initialiseDisplacement();

            //- Update global patches
            virtual void update(const bool force = false);

            //- Apply topological changes
            virtual bool changeMesh();

            //- Initialise the mesh
            virtual void initialiseMesh(const IterType iter);

            //- Initialise the region fields
            virtual void initialiseFields();

            //- Initialise the solver
            virtual void initialise();

            //- Move the mesh
            virtual bool moveMesh(const IterType iter);

            //- Solve the model
            virtual void solve();

            //- Clear the model
            virtual void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
