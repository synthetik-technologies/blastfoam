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
    Foam::immersedBoundaryObjectRestraints::axialAngularSpring

Description
    immersedBoundaryObjectRestraints model.  Axial angular spring with moment
    values drawn from an interpolation table.  Linear damping.

SourceFiles
    axialAngularSpring.C

\*---------------------------------------------------------------------------*/

#ifndef axialAngularSpring_H
#define axialAngularSpring_H

#include "immersedBoundaryObjectRestraint.H"
#include "point.H"
#include "tensor.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace immersedBoundaryObjectRestraints
{

/*---------------------------------------------------------------------------*\
                  Class axialAngularSpring Declaration
\*---------------------------------------------------------------------------*/

class axialAngularSpring
:
    public immersedBoundaryObjectRestraint
{
    // Private Data

        //- Reference orientation where there is no moment
        tensor refQ_;

        //- Global unit axis around which the motion is sprung
        vector axis_;

        //- Spring moment interpolation table, depending on angleFormat
        autoPtr<Function1<scalar>> moment_;

        //- Boolean stating whether the angle around the axis needs to
        //  be converted to degrees before asking the
        //  interpolationTable for a moment value
        bool convertToDegrees_;

        //- Damping coefficient (Nms/rad)
        scalar damping_;


public:

    //- Runtime type information
    TypeName("axialAngularSpring");


    // Constructors

        //- Construct from components
        axialAngularSpring
        (
            const polyMesh& mesh,
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<immersedBoundaryObjectRestraint> clone() const
        {
            return autoPtr<immersedBoundaryObjectRestraint>
            (
                new axialAngularSpring(*this)
            );
        }


    //- Destructor
    virtual ~axialAngularSpring();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const immersedBoundaryObject& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
