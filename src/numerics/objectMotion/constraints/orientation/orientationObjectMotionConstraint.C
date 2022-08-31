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

\*---------------------------------------------------------------------------*/

#include "orientationObjectMotionConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "movingObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionConstraints
{
    defineTypeNameAndDebug(orientation, 0);

    addToRunTimeSelectionTable
    (
        objectMotionConstraint,
        orientation,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::orientation::orientation
(
    const word& name,
    const dictionary& dict,
    const movingObject& motion
)
:
    objectMotionConstraint(name, dict, motion)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::orientation::~orientation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionConstraints::orientation::
constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::objectMotionConstraints::orientation::
constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(3, Zero)));
}


bool Foam::objectMotionConstraints::orientation::read
(
    const dictionary& dict
)
{
    objectMotionConstraint::read(dict);

    return true;
}


void Foam::objectMotionConstraints::orientation::write
(
    Ostream& os
) const
{
}

// ************************************************************************* //
