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
    Foam::movingObject

Description
    Base class for a moving object


SourceFiles
    movingObject.C

\*---------------------------------------------------------------------------*/

#ifndef movingObject_H
#define movingObject_H

#include "point.H"
#include "vector.H"
#include "tensor.H"
#include "diagTensor.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class movingObject Declaration
\*---------------------------------------------------------------------------*/

class movingObject
{
public:

    TypeName("movingObject");

    // Constructors

        //- Null constructor
        movingObject();

    //- Destructor
    virtual ~movingObject();


    // Member Functions

        // Access

            //- Print information
            virtual bool report() const = 0;

            // Difference between center of mass and centre of rotation
            virtual vector momentArm() const = 0;

            //- Return the current centre of rotation
            virtual point centreOfRotation() const = 0;

            //- Return the initial centre of mass
            virtual point initialCentreOfMass() const = 0;

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

            //- Return the current velocity at a point
            virtual vector v(const point& pt) const = 0;

            //- Return the inertia tensor
            virtual diagTensor momentOfInertia() const = 0;

            //- Transform the given initial state point by the current motion
            //  state
            virtual point transform(const point&) const = 0;

            //- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
            //  and return the rotated Q and pi as a tuple
            virtual Tuple2<tensor, vector> rotate
            (
                const tensor& Q0,
                const vector& pi,
                const scalar deltaT
            ) const = 0;


        // Update state

            //- Update and relax accelerations from the force and torque
            virtual void updateAcceleration
            (
                const vector& fGlobal,
                const vector& tauGlobal
            ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
