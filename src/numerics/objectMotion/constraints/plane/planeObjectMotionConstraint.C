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

#include "planeObjectMotionConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "movingObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionConstraints
{
    defineTypeNameAndDebug(plane, 0);

    addToRunTimeSelectionTable
    (
        objectMotionConstraint,
        plane,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::plane::plane
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

Foam::objectMotionConstraints::plane::~plane()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionConstraints::plane::setCentreOfRotation
(
    point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::objectMotionConstraints::plane::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.applyConstraint(normal_);
}


void Foam::objectMotionConstraints::plane::constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, normal_)));
}


bool Foam::objectMotionConstraints::plane::read
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

    coeffDict_.lookup("normal") >> normal_;

    return true;
}


void Foam::objectMotionConstraints::plane::write
(
    Ostream& os
) const
{
    writeEntry(os, "centreOfRotation", centreOfRotation_);
    writeEntry(os, "normal", normal_);
}

// ************************************************************************* //
