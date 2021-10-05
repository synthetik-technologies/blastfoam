/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::dynamicMotionSolverListBlastFvMesh

Description
    Foam::dynamicMotionSolverListBlastFvMesh

SourceFiles
    dynamicMotionSolverListBlastFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverListBlastFvMesh_H
#define dynamicMotionSolverListBlastFvMesh_H

#include "dynamicBlastFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
            Class dynamicMotionSolverListBlastFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverListBlastFvMesh
:
    public dynamicBlastFvMesh
{
    // Private data

        PtrList<motionSolver> motionSolvers_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicMotionSolverListBlastFvMesh(const dynamicMotionSolverListBlastFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverListBlastFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverListFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverListBlastFvMesh(const IOobject& io);


    //- Destructor
    ~dynamicMotionSolverListBlastFvMesh();


    // Member Functions

        //- Dummy update function which does not change the mesh
        virtual bool update();

        //- Refine the mesh
        virtual bool refine(const bool correctError = true);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
