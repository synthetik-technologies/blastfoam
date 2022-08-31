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

\*---------------------------------------------------------------------------*/

#include "axisObjectMotionConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionConstraints
{
    defineTypeNameAndDebug(axis, 0);

    addToRunTimeSelectionTable
    (
        objectMotionConstraint,
        axis,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::axis::axis
(
    const word& name,
    const dictionary& dict,
    const movingObject& motion
)
:
    objectMotionConstraint(name, dict, motion),
    axis_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::axis::~axis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionConstraints::axis::constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::objectMotionConstraints::axis::constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
}


bool Foam::objectMotionConstraints::axis::read
(
    const dictionary& dict
)
{
    objectMotionConstraint::read(dict);

    coeffDict_.lookup("axis") >> axis_;

    scalar magFixedAxis(mag(axis_));

    if (magFixedAxis > vSmall)
    {
        axis_ /= magFixedAxis;
    }
    else
    {
        FatalErrorInFunction
            << "axis has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::objectMotionConstraints::axis::write
(
    Ostream& os
) const
{
    writeEntry(os, "axis", axis_);
}

// ************************************************************************* //
