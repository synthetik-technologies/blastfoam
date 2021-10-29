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
    Foam::MovingImmersedBoundaryObject

Description


SourceFiles
    MovingImmersedBoundaryObjectI.H
    MovingImmersedBoundaryObject.C

\*---------------------------------------------------------------------------*/

#ifndef MovingImmersedBoundaryObject_H
#define MovingImmersedBoundaryObject_H

#include "movingImmersedBoundaryObjectState.H"
#include "pointField.H"
#include "immersedBoundaryObject.H"
#include "immersedBoundaryObjectRestraint.H"
#include "immersedBoundaryObjectConstraint.H"
#include "Tuple2.H"
#include "dictionaryEntry.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class immersedBoundaryObjectMotionSolver;

/*---------------------------------------------------------------------------*\
                      Class MovingImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class MovingImmersedBoundaryObject
:
    public ImmersedType
{
    // Private Data

        //- Motion state data object
        movingImmersedBoundaryObjectState motionState_;

        //- Motion state data object for previous time-step
        movingImmersedBoundaryObjectState motionState0_;

        //- Motion restraints
        PtrList<immersedBoundaryObjectRestraint> restraints_;

        //- Motion constraints
        PtrList<immersedBoundaryObjectConstraint> constraints_;

        //- Translational constraint tensor
        tensor tConstraints_;

        //- Rotational constraint tensor
        tensor rConstraints_;

        //- Orientation of centre of rotation
        vector initialCentreOfRotation_;

        //- Moment of inertia of the body in reference configuration
        //  (Q = I)
        diagTensor momentOfInertia_;

        //- Acceleration relaxation coefficient
        scalar aRelax_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;

        //- Motion solver
        autoPtr<immersedBoundaryObjectMotionSolver> solver_;

        //- Initial points
        pointField points0_;


    // Private Member Functions

        //- Calculate the rotation tensor around the body reference
        //  frame x-axis by the given angle
        inline tensor rotationTensorX(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame y-axis by the given angle
        inline tensor rotationTensorY(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame z-axis by the given angle
        inline tensor rotationTensorZ(scalar deltaT) const;

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


        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return the restraints
            inline const PtrList<immersedBoundaryObjectRestraint>&
                restraints() const;

            //- Return the constraints
            inline const PtrList<immersedBoundaryObjectConstraint>&
                constraints() const;

            //- Return the orientation
            inline const tensor& Q() const;

            //- Return the current acceleration
            inline const vector& a() const;

            //- Return the current angular momentum
            inline const vector& pi() const;

            //- Return the current torque
            inline const vector& tau() const;


        // Edit

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to vector
            inline vector& v();

            //- Return non-const access to acceleration
            inline vector& a();

            //- Return non-const access to angular momentum
            inline vector& pi();

            //- Return non-const access to torque
            inline vector& tau();


public:

    // Constructors

        //- Construct from constant and state dictionaries
        MovingImmersedBoundaryObject
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~MovingImmersedBoundaryObject();


    // Member Functions

        // Access

            //- Finish setting up the object
            virtual void initialize();

            //- Is the object moving
            virtual bool moving() const
            {
                return true;
            }

            //- Return the motion state
            inline const movingImmersedBoundaryObjectState& state() const;

            virtual vector momentArm() const;

            //- Return the initial points
            inline const pointField& points0() const;

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

            //- Return the current velocity
            virtual vector v() const;

            //- Return the inertia tensor
            virtual diagTensor momentOfInertia() const;


        // Edit

            //- Store the motion state at the beginning of the time-step
            inline void newTime();

            //- Return non-const access to the centre of rotation
            inline point& centreOfRotation();


        // Constraints and Restraints

            //- Add restraints to the motion, public to allow external
            //  addition of restraints after construction
            void addRestraints(const dictionary& dict);

            //- Add restraints to the motion, public to allow external
            //  addition of restraints after construction
            void addConstraints(const dictionary& dict);


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
            );

            //- Update
            virtual void update();

            //- Solve
            virtual void solve();

            //- Report the status of the motion
            virtual void status() const;


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
    #include "MovingImmersedBoundaryObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "MovingImmersedBoundaryObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
