/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025 Synthetik Applied Technologies
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
    for more detonationails.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::errorTypes::field

Description
    Computes detonation velocity based on the time when a target pressure is
    reached at two different locations

SourceFiles
    detonationVelocityError.C

\*---------------------------------------------------------------------------*/

#ifndef detonationVelocityError_H
#define detonationVelocityError_H

#include "errorType.H"
#include "runTimeSelectionTables.H"


namespace Foam
{
namespace errorTypes
{

/*---------------------------------------------------------------------------*\
                            Class detonationVelocity Declaration
\*---------------------------------------------------------------------------*/

class detonationVelocity
:
    public errorType
{
protected:
    // Protected data

        //- Name of the field
        const word pName_;

        //- Probe locations
        vector point1_;
        vector point2_;

        //- Saved detonation velocity
        scalar vDet_;

        //- Time when first probe is reached
        scalar startTime_;

        //- Time when second probe is reached
        scalar endTime_;

        //- Stop the simulation when the second probe is reached
        bool stopWhenReached_;

        //- Pressure to mark "reached"
        scalar pTarget_;


public:

    //- Runtime type information
    TypeName("detonationonationVelocity");

    // Constructor

        // Construct from dictionary
        detonationVelocity
        (
            const Time& runTime,
            const dictionary& dict,
            const word& region
        );

    //- Clone the object
    virtual autoPtr<errorType> clone() const
    {
        return autoPtr<errorType>(new detonationVelocity(*this));
    }


    //- Destructor
    virtual ~detonationVelocity();


    //Member functions

        //- Return the field name
        virtual const word& fieldName() const
        {
            return pName_;
        }

        //- Update the error
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
