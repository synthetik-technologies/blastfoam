/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020   Jeff Heylmun:  Set old cell volumes after refinement is done.
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
    Foam::movingAdaptiveFvMesh

SourceFiles
    movingAdaptiveFvMesh.C

Description
    Mesh movement is solved with motion solvers, and adaptive mesh refinement
    is also included.

Additional files
    adaptiveFvMesh.H

\*---------------------------------------------------------------------------*/

#ifndef movingAdaptiveFvMesh_H
#define movingAdaptiveFvMesh_H

#include "adaptiveFvMesh.H"
#include "hexRef.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class movingAdaptiveFvMesh Declaration
\*---------------------------------------------------------------------------*/

class motionSolver;

class movingAdaptiveFvMesh
:
    public adaptiveFvMesh
{
    // Private data

        autoPtr<motionSolver> motionPtr_;

        //- Optional list of vectorFields to update for mesh motion
        //  For modern solvers using Uf and correctPhi to update the flux
        //  after motion it is not necessary to specify a "velocityFields" list
        velocityMotionCorrection velocityMotionCorrection_;


    // Private Member Functions

        //- Overload update mesh to include other methods
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Disallow default bitwise copy construct
        movingAdaptiveFvMesh(const movingAdaptiveFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const movingAdaptiveFvMesh&);



public:


     //- Runtime type information
    TypeName("movingAdaptiveFvMesh");


    // Constructors

        //- Construct from IOobject
        movingAdaptiveFvMesh(const IOobject& io);


    //- Destructor
    ~movingAdaptiveFvMesh();


    // Member Functions

        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


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
