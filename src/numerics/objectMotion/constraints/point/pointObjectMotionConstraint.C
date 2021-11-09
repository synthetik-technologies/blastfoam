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

#include "pointObjectMotionConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "movingObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionConstraints
{
    defineTypeNameAndDebug(point, 0);

    addToRunTimeSelectionTable
    (
        objectMotionConstraint,
        point,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::point::point
(
    const word& name,
    const dictionary& dict,
    const movingObject& motion
)
:
    objectMotionConstraint(name, dict, motion),
    centreOfRotation_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionConstraints::point::~point()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionConstraints::point::setCentreOfRotation
(
    Foam::point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::objectMotionConstraints::point::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(3, Zero)));
}


void Foam::objectMotionConstraints::point::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::objectMotionConstraints::point::read
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

    return true;
}


void Foam::objectMotionConstraints::point::write
(
    Ostream& os
) const
{
    writeEntry(os, "centreOfRotation", centreOfRotation_);
}

// ************************************************************************* //
