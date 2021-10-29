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
    Foam::TemperatureIndependentImmersedBoundaryObject

Description


SourceFiles
    TemperatureIndependentImmersedBoundaryObject.C

\*---------------------------------------------------------------------------*/

#ifndef TemperatureIndependentImmersedBoundaryObject_H
#define TemperatureIndependentImmersedBoundaryObject_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class TemperatureIndependentImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class TemperatureIndependentImmersedBoundaryObject
:
    public ImmersedType
{

public:

    // Constructors

        //- Construct from constant and state dictionaries
        TemperatureIndependentImmersedBoundaryObject
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~TemperatureIndependentImmersedBoundaryObject();


    // Member Functions

        //- Is the object dependent on temperature
        virtual bool temperatureDependent() const
        {
            return false;
        }

        // Update state

            //- Report the status of the motion
            virtual void status() const;


        // Transformations

            //- Move boundary points
            virtual void movePoints();


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
    #include "TemperatureIndependentImmersedBoundaryObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
