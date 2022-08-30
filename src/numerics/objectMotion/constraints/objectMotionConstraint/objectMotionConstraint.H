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

Namespace
    Foam::objectMotionConstraints

Description
    Namespace for six DoF motion constraints


Class
    Foam::objectMotionConstraint

Description
    Base class for defining constraints for sixDoF motions

SourceFiles
    objectMotionConstraint.C
    objectMotionConstraintNew.C

\*---------------------------------------------------------------------------*/

#ifndef objectMotionConstraint_H
#define objectMotionConstraint_H

#include "dictionary.H"
#include "autoPtr.H"
#include "point.H"
#include "pointConstraint.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class movingObject;

/*---------------------------------------------------------------------------*\
                Class objectMotionConstraint Declaration
\*---------------------------------------------------------------------------*/

class objectMotionConstraint
{

protected:

    // Protected data

        //- Name of the constraint
        word name_;

        //- Constraint model specific coefficient dictionary
        dictionary coeffDict_;

        //- Moving object
        const movingObject& motion_;


public:

    //- Runtime type information
    TypeName("objectMotionConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectMotionConstraint,
            motion,
            (
                const word& name,
                const dictionary& dict,
                const movingObject& motion
            ),
            (name, dict, motion)
        );


    // Constructors

        //- Construct from the sDoFRBMCDict dictionary and Time
        objectMotionConstraint
        (
            const word& name,
            const dictionary& dict,
            const movingObject& motion
        );

        //- Construct and return a clone
        virtual autoPtr<objectMotionConstraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMCDict dictionary and Time
        static autoPtr<objectMotionConstraint> New
        (
            const word& name,
            const dictionary& dict,
            const movingObject& motion
        );


    //- Destructor
    virtual ~objectMotionConstraint();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        //- Set the centre of rotation if not the centre of mass
        virtual void setCentreOfRotation(point&) const
        {}

        //- Apply and accumulate translational constraints
        virtual void constrainTranslation(pointConstraint&) const = 0;

        //- Apply and accumulate rotational constraints
        virtual void constrainRotation(pointConstraint&) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        // Access

            // Return access to sDoFRBMCCoeffs
            inline const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
