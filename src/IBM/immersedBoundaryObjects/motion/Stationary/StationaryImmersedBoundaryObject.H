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
public:

    // Constructors

        //- Construct from constant and state dictionaries
        StationaryImmersedBoundaryObject
        (
            const polyPatch& patch,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~StationaryImmersedBoundaryObject();


    // Member Functions

        // Access

            //- Is the object moving
            virtual bool moving() const
            {
                return false;
            }

            //- Return the centre
            virtual point centre() const;

            //- Return the mass of the object
            virtual scalar mass() const;

            //- Return the rotated inertia tensor
            virtual diagTensor momentOfInertia() const;


        // Update state

            //- Report the status of the motion
            virtual void status(const bool print = true) const;


        // Transformations

            //- Return the velocity at a position
            virtual vector v(const vector& pt) const;

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity(const pointField& pt) const;

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity() const;

            //- Return the orientation tensor, Q.
            //  globalVector = Q & bodyLocalVector
            //  bodyLocalVector = Q.T() & globalVector
            virtual tensor orientation() const;

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
