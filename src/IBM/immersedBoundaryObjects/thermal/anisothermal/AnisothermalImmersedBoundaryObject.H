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
    Foam::AnisothermalImmersedBoundaryObject

Description


SourceFiles
    AnisothermalImmersedBoundaryObject.C

\*---------------------------------------------------------------------------*/

#ifndef AnisothermalImmersedBoundaryObject_H
#define AnisothermalImmersedBoundaryObject_H

#include "solidBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class AnisothermalImmersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class AnisothermalImmersedBoundaryObject
:
    public ImmersedType
{
    // Private Data

        //- Solid thermodynamic model
        autoPtr<solidBlastThermo> thermo_;



public:

    // Constructors

        //- Construct from constant and state dictionaries
        AnisothermalImmersedBoundaryObject
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~AnisothermalImmersedBoundaryObject();


    // Member Functions

        //- Is the object dependent on temperature
        virtual bool temperatureDependent() const
        {
            return true;
        }

        //- Initialize
        virtual void initialize();

        // Update state

            //- Update
            virtual void solve();

            //- Report the status of the motion
            virtual void status(const bool print = true) const;


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
    #include "AnisothermalImmersedBoundaryObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
