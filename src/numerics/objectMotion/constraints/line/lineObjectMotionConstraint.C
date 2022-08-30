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

#include "lineObjectMotionConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "movingObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionConstraints
{
    defineTypeNameAndDebug(line, 0);

    addToRunTimeSelectionTable
    (
        objectMotionConstraint,
        line,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::line::line
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

Foam::objectMotionConstraints::line::~line()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionConstraints::line::setCentreOfRotation
(
    point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::objectMotionConstraints::line::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, direction_)));
}


void Foam::objectMotionConstraints::line::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::objectMotionConstraints::line::read
(
    const dictionary& dict
)
{
    objectMotionConstraint::read(dict);

    centreOfRotation_ = coeffDict_.lookupOrDefault
    (
        "centreOfRotation",
        motion_.initialCentreOfMass()
    );

    coeffDict_.lookup("direction") >> direction_;

    scalar magDir(mag(direction_));

    if (magDir > vSmall)
    {
        direction_ /= magDir;
    }
    else
    {
        FatalErrorInFunction
            << "line direction has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::objectMotionConstraints::line::write
(
    Ostream& os
) const
{
    writeEntry(os, "centreOfRotation", centreOfRotation_);
    writeEntry(os, "direction", direction_);
}

// ************************************************************************* //
