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

#include "objectMotionState.H"
#include "objectMotionSolver.H"
#include "pointField.H"
#include "immersedBoundaryObject.H"
#include "movingObject.H"
#include "objectMotionState.H"
#include "objectMotionRestraint.H"
#include "objectMotionConstraint.H"
#include "Tuple2.H"
#include "dictionaryEntry.H"
#include "PtrListDictionary.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class objectMotionSolver;

/*---------------------------------------------------------------------------*\
                      Class MovingImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class MovingImmersedBoundaryObject
:
    public movingObject,
    public ImmersedType
{
    // Private Data

        //- Motion state data object
        objectMotionState motionState_;

        //- Motion state data object for previous time-step
        objectMotionState motionState0_;

        //- Motion restraints
        PtrList<objectMotionRestraint> restraints_;

        //- Motion constraints
        PtrList<objectMotionConstraint> constraints_;

        //- Translational constraint tensor
        tensor tConstraints_;

        //- Rotational constraint tensor
        tensor rConstraints_;

        //- Orientation of centre of rotation
        vector initialCentreOfRotation_;

        //- Orientation of centre of mass
        vector initialCentreOfMass_;

        //- Moment of inertia of the body in reference configuration
        //  (Q = I)
        mutable diagTensor momentOfInertia_;

        //- Acceleration relaxation coefficient
        scalar aRelax_;

        //- Motion solver
        autoPtr<objectMotionSolver> solver_;

        //- Initial points
        pointField points0_;

        //- Mass of object
        scalar mass_;

        //- Gravitational acceleration
        vector g_;

        //- Scale factor for the mesh
        scalar scale_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;

        //- Does this object need to be initialized
        bool needToInitialize_;

        //- Optional log file
        mutable autoPtr<OFstream> outputPtr_;

    // Private Member Functions

        //- Initialize the variables based on the mesh shape
        void initializeState();

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

        //- Return the scale factor for the mesh
        scalar scale() const;

        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return the restraints
            inline const PtrList<objectMotionRestraint>& restraints() const;

            //- Return the constraints
            inline const PtrList<objectMotionConstraint>& constraints() const;

            //- Return the orientation
            inline const tensor& Q() const;

            //- Return the current acceleration
            inline const vector& a() const;

            //- Return the current angular momentum
            inline const vector& pi() const;

            //- Return the current torque
            inline const vector& tau() const;

            //- Return the report Switch
            inline bool report() const;


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
            const polyPatch& patch,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~MovingImmersedBoundaryObject();


    // Member Functions

        // Access

            //- Is the object moving
            virtual bool moving() const
            {
                return true;
            }

            //- Return the mass of the object
            virtual scalar mass() const;

            //- Return the motion state
            inline const objectMotionState& state() const;

            virtual vector momentArm() const;

            //- Return the initial points
            inline const pointField& points0() const;

            //- Return the initial centre of rotation
            virtual point initialCentreOfRotation() const;

            //- Return the current centre of rotation
            virtual point centreOfRotation() const;

            //- Return the initial centre of mass
            virtual point initialCentreOfMass() const;

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

            //- Return the current velocity at a point
            virtual vector v(const point& pt) const;

            //- Return the current velocity
            virtual vector v() const;

            //- Return the inertia tensor
            virtual diagTensor momentOfInertia() const;


        // Edit

            //- Store the motion state at the beginning of the time-step
            inline void newTime();

            //- Return non-const access to the centre of rotation
            inline point& centreOfRotation();

            //- Return non-const access to the centre of rotation
            virtual point centre() const
            {
                return centreOfRotation();
            }


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
            virtual void status(const bool print = true) const;


        // Transformations

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity(const pointField& pt) const;

            //- Return the velocity at the faceCentres
            virtual tmp<vectorField> velocity() const;

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

        virtual void logHeader(Ostream&) const;
        virtual void log(Ostream&) const;

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
