/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::motionImmersedBoundaryObject

Description


SourceFiles
    motionImmersedBoundaryObjectI.H
    motionImmersedBoundaryObject.C

\*---------------------------------------------------------------------------*/

#ifndef motionImmersedBoundaryObject_H
#define motionImmersedBoundaryObject_H

#include "pointField.H"
#include "immersedBoundaryObject.H"
#include "Tuple2.H"
#include "dictionaryEntry.H"
#include "PtrListDictionary.H"
#include "immersedBoundaryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class immersedBoundaryObjectMotionSolver;

/*---------------------------------------------------------------------------*\
                      Class motionImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/
\
class motionImmersedBoundaryObject
{



public:

    // Constructors

        //- Construct from constant and state dictionaries
        motionImmersedBoundaryObject
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~motionImmersedBoundaryObject();


    // Member Functions

        // Access

            //- Is the object motion
            virtual bool motion() const = 0;

            // Difference between center of mass and centre of rotation
            virtual vector momentArm() const = 0;

            //- Return the current centre of rotation
            virtual point centreOfRotation() const = 0;

            //- Return the current centre of mass
            virtual point centreOfMass() const = 0;

            //- Return the translational constraint tensor
            virtual tensor tConstraints() const = 0;

            //- Return the rotational constraint tensor
            virtual tensor rConstraints() const = 0;

            //- Return the orientation tensor, Q.
            //  globalVector = Q & bodyLocalVector
            //  bodyLocalVector = Q.T() & globalVector
            virtual tensor orientation() const = 0;

            //- Return the angular velocity in the global frame
            virtual vector omega() const = 0;

            //- Return the current velocity
            virtual vector v() const = 0;

            //- Return the inertia tensor
            virtual diagTensor momentOfInertia() const = 0;


        // Update state

            //- Symplectic integration of velocities, orientation and position.
            //  Changes to Crank-Nicolson integration for subsequent iterations.
            virtual bool update
            (
                bool firstIter,
                const vector& fGlobal,
                const vector& tauGlobal,
                scalar deltaT,
                scalar deltaT0
            ) = 0;

            //- Report the status of the motion
            virtual void status() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
