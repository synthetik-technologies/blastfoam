/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    frictionLaw

Description
    Friction law e.g. Coulomb, Orowan, Shaw, etc.
    Given the pressure (or other variable for more complex laws),
    the friction law returns the shear traction.

SourceFiles
    frictionLaw.C
    newFrictionLaw.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef frictionLaw_H
#define frictionLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class frictionContactModel;

/*---------------------------------------------------------------------------*\
                         Class frictionLaw Declaration
\*---------------------------------------------------------------------------*/

class frictionLaw
{
    // Private data

        //- Reference to the frictionContactModel
        const frictionContactModel& frictionModel_;

        //- Friction law dictionary
        dictionary frictionLawDict_;


    // Private Member Functions

        //- Disallow copy construct
        //frictionLaw(const frictionLaw&);

        //- Disallow default bitwise assignment
        void operator=(const frictionLaw&);


public:

    //- Runtime type information
    TypeName("frictionLaw");

    // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            frictionLaw,
            dictionary,
            (
                const word name,
                const frictionContactModel& fricModel,
                const dictionary& dict
            ),
            (name, fricModel, dict)
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<frictionLaw> New
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        frictionLaw
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        frictionLaw(const frictionLaw& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const = 0;


    // Destructor

        virtual ~frictionLaw()
        {}


    // Member Functions

        //- Return reference to friction contact model
        const frictionContactModel& frictionModel() const
        {
            return frictionModel_;
        }

        //- Return friction law dict
        const dictionary& frictionLawDict() const
        {
            return frictionLawDict_;
        }

        //- Return traction that would cause slipping
        virtual scalar slipTraction(const scalar pressure);

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
            const scalar pressure, const vector& slipDir
        );

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
                const scalar contactPressure,         // Contact pressure
                const vector& faceSlip,               // Slip vector
                const vector& slaveFaceVelocity,      // Velocity of slave face
                const vector& masterFaceVelocity      // Velocity of master face
        );

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
            const scalar contactPressure,         // Contact pressure
            const vector& faceSlip,               // Slip vector
            const vector& slaveFaceVelocity,      // Velocity of slave face
            const vector& masterFaceVelocity,     // Velocity of master face
            const label slavePatchIndex,          // Slave patch index
            const label faceIndex                 // Local slave face ID
        );

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
