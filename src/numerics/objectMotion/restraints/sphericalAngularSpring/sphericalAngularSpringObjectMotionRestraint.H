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
    Foam::objectMotionRestraints::sphericalAngularSpring

Description
    objectMotionRestraints model.  Spherical angular spring.

SourceFiles
    sphericalAngularSpringObjectMotionRestraint.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalAngularSpringObjectMotionRestraint_H
#define sphericalAngularSpringObjectMotionRestraint_H

#include "objectMotionRestraint.H"
#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectMotionRestraints
{

/*---------------------------------------------------------------------------*\
                   Class sphericalAngularSpring Declaration
\*---------------------------------------------------------------------------*/

class sphericalAngularSpring
:
    public objectMotionRestraint
{
    // Private Data

        //- Reference orientation where there is no moment
        tensor refQ_;

        //- Spring stiffness coefficient (Nm/rad)
        scalar stiffness_;

        //- Damping coefficient (Nms/rad)
        scalar damping_;


public:

    //- Runtime type information
    TypeName("sphericalAngularSpring");


    // Constructors

        //- Construct from components
        sphericalAngularSpring
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<objectMotionRestraint> clone() const
        {
            return autoPtr<objectMotionRestraint>
            (
                new sphericalAngularSpring(*this)
            );
        }


    //- Destructor
    virtual ~sphericalAngularSpring();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const movingObject& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
