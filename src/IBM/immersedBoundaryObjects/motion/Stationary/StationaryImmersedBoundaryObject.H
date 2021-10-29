/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
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
    Foam::StationaryImmersedBoundaryObject

Description


SourceFiles
    StationaryImmersedBoundaryObject.C

\*---------------------------------------------------------------------------*/

#ifndef StationaryImmersedBoundaryObject_H
#define StationaryImmersedBoundaryObject_H

#include "pointField.H"
#include "immersedBoundaryObject.H"
#include "dictionaryEntry.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class immersedBoundaryObjectMotionSolver;

/*---------------------------------------------------------------------------*\
                      Class StationaryImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class StationaryImmersedBoundaryObject
:
    public ImmersedType
{
    // Private data

        //- Orientation of the object
        tensor orientation_;

    // Private Functions

        //- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
        //  and return the rotated Q and pi as a tuple
        virtual Tuple2<tensor, vector> rotate
        (
            const tensor& Q0,
            const vector& pi,
            const scalar deltaT
        ) const;

        //- Apply the restraints to the object
        virtual void applyRestraints();

        //- Update and relax accelerations from the force and torque
        virtual void updateAcceleration
        (
            const vector& fGlobal,
            const vector& tauGlobal
        );


public:

    // Constructors

        //- Construct from constant and state dictionaries
        StationaryImmersedBoundaryObject
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~StationaryImmersedBoundaryObject();


    // Member Functions

        // Access

            //- Finish setting up the object
            virtual void initialize();

            //- Is the object moving
            virtual bool moving() const
            {
                return false;
            }

            //- Return the current centre of rotation
            virtual point centreOfRotation() const;

            //- Return the current centre of mass
            virtual point centreOfMass() const;

            //- Return the translational constraint tensor
            virtual tensor tConstraints() const;

            //- Return the rotational constraint tensor
            virtual tensor rConstraints() const;

            //- Return the orientation tensor, Q.
            //  globalVector = Q & bodyLocalVector
            //  bodyLocalVector = Q.T() & globalVector
            virtual tensor orientation() const;

            //- Return the angular velocity in the global frame
            virtual vector omega() const;

            virtual vector momentArm() const;

            //- Return the current velocity
            virtual vector v() const;

            //- Return the inertia tensor
            virtual diagTensor momentOfInertia() const;


        // Update state

            //- Report the status of the motion
            virtual void status() const;

            //- Store the motion state at the beginning of the time-step
            virtual void newTime()
            {}


        // Transformations

            //- Return the velocity at a position
            virtual vector velocity(const vector& pt) const;

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity(const pointField& pt) const;

            //- Transform the given initial state point by the current motion
            //  state
            virtual point transform(const point&) const;

            //- Transform the given initial state pointField by the current
            //  motion state
            virtual tmp<pointField> transform(const pointField&) const;

            //- Transform the given initial state pointField by the current
            //  motion state
            virtual point inverseTransform(const point&) const;

            //- Transform the given initial state pointField by the current
            //  motion state
            virtual tmp<pointField> inverseTransform(const pointField&) const;

            //- Move boundary points
            virtual void movePoints();


        //- Write
        virtual void write(Ostream&) const;
        virtual void write(dictionary&) const;

        //- Read coefficients dictionary and update system parameters,
        //  constraints and restraints but not the current state
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "StationaryImmersedBoundaryObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
