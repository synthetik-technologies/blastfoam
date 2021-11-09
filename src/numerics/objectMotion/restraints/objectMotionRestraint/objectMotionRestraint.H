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
    Foam::objectMotionRestraints

Description
    Namespace for six DoF motion restraints


Class
    Foam::objectMotionRestraint

Description
    Base class for defining restraints for sixDoF motions

SourceFiles
    objectMotionRestraint.C
    objectMotionRestraintNew.C

\*---------------------------------------------------------------------------*/

#ifndef objectMotionRestraint_H
#define objectMotionRestraint_H

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "polyMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class movingObject;

/*---------------------------------------------------------------------------*\
                Class objectMotionRestraint Declaration
\*---------------------------------------------------------------------------*/

class objectMotionRestraint
{

protected:

    // Protected data

        //- Name of the restraint
        word name_;

        //- Restraint model specific coefficient dictionary
        dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("objectMotionRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectMotionRestraint,
            motion,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );


    // Constructors

        //- Construct from a dictionary and mesh
        objectMotionRestraint
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<objectMotionRestraint> clone() const = 0;


    // Selectors

        //- Select constructed from a dictionary and Time
        static autoPtr<objectMotionRestraint> New
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~objectMotionRestraint();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const movingObject& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        // Access

            // Return access to sDoFRBMRCoeffs
            inline const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

        //- Write
        virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
