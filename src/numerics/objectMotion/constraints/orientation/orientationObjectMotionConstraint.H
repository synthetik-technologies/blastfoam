/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::objectMotionConstraints::orientation

Description
    Orientation constraint:
        fixed in global space.

SourceFiles
    orientationObjectMotionConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef orientationObjectMotionConstraint_H
#define orientationObjectMotionConstraint_H

#include "objectMotionConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectMotionConstraints
{

/*---------------------------------------------------------------------------*\
                       Class orientation Declaration
\*---------------------------------------------------------------------------*/

class orientation
:
    public objectMotionConstraint
{

public:

    //- Runtime type information
    TypeName("orientation");


    // Constructors

        //- Construct from components
        orientation
        (
            const word& name,
            const dictionary& dict,
            const movingObject& motion
        );

        //- Construct and return a clone
        virtual autoPtr<objectMotionConstraint> clone() const
        {
            return autoPtr<objectMotionConstraint>
            (
                new orientation(*this)
            );
        }


    //- Destructor
    virtual ~orientation();


    // Member Functions

        //- Apply and accumulate translational constraints
        virtual void constrainTranslation(pointConstraint&) const;

        //- Apply and accumulate rotational constraints
        virtual void constrainRotation(pointConstraint&) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectMotionConstraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
