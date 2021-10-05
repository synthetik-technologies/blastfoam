/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::dynamicMotionSolverBlastFvMesh

Description
    The dynamicMotionSolverBlastFvMesh

SourceFiles
    dynamicMotionSolverBlastFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverBlastFvMesh_H
#define dynamicMotionSolverBlastFvMesh_H

#include "dynamicBlastFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                  Class dynamicMotionSolverBlastFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverBlastFvMesh
:
    public dynamicBlastFvMesh
{
    // Private data

        autoPtr<motionSolver> motionPtr_;

        //- Optional list of vectorFields to update for mesh motion
        //  For modern solvers using Uf and correctPhi to update the flux
        //  after motion it is not necessary to specify a "velocityFields" list
        velocityMotionCorrection velocityMotionCorrection_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicMotionSolverBlastFvMesh(const dynamicMotionSolverBlastFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverBlastFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverBlastFvMesh(const IOobject& io);


    //- Destructor
    ~dynamicMotionSolverBlastFvMesh();


    // Member Functions

        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Refine the mesh
        virtual bool refine(const bool correctError = true);


        // Write

            //- Write the underlying polyMesh and other data
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
