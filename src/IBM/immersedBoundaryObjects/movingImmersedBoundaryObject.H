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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class motionImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/

class motionImmersedBoundaryObject
{
public:

    TypeName("motionImmersedBoundaryObject");

    // Constructors

        //- Construct from constant and state dictionaries
        motionImmersedBoundaryObject();

    //- Destructor
    virtual ~motionImmersedBoundaryObject();


    // Member Functions

        // Access

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

            //- Return the velocity at a list of position
            virtual vector v(const point& pt) const = 0;

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity(const pointField& pt) const = 0;

            //- Return the current velocity
            virtual tmp<vectorField> velocity() const = 0;

            //- Return the inertia tensor
            virtual diagTensor momentOfInertia() const = 0;

            //- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
            //  and return the rotated Q and pi as a tuple
            virtual Tuple2<tensor, vector> rotate
            (
                const tensor& Q0,
                const vector& pi,
                const scalar deltaT
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
