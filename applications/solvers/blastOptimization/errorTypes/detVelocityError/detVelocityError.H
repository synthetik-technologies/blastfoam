/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
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
    Foam::errorTypes::field

Description
    Compute error in field estimates. Can use operations.

SourceFiles
    detVelocityError.C

\*---------------------------------------------------------------------------*/

#ifndef detVelocityError_H
#define detVelocityError_H

#include "errorType.H"
#include "runTimeSelectionTables.H"


namespace Foam
{
namespace errorTypes
{

/*---------------------------------------------------------------------------*\
                            Class detVelocity Declaration
\*---------------------------------------------------------------------------*/

class detVelocity
:
    public errorType
{
protected:
    // Protected data

        //- Name of the field
        const word pName_;

        vector point1_;
        vector point2_;

        scalar vDet_;

        scalar startTime_;
        scalar endTime_;
        bool stopWhenReached_;

        scalar pTarget_;


public:

    //- Runtime type information
    TypeName("detonationVelocity");

    // Constructor

        // Construct from dictionary
        detVelocity
        (
            const Time& runTime,
            const dictionary& dict,
            const word& region
        );

    //- Clone the object
    virtual autoPtr<errorType> clone() const
    {
        return autoPtr<errorType>(new detVelocity(*this));
    }


    //- Destructor
    virtual ~detVelocity();


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
