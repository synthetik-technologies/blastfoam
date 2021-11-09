/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "sphericalAngularDamperObjectMotionRestraint.H"
#include "addToRunTimeSelectionTable.H"
#include "movingObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionRestraints
{
    defineTypeNameAndDebug(sphericalAngularDamper, 0);

    addToRunTimeSelectionTable
    (
        objectMotionRestraint,
        sphericalAngularDamper,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionRestraints::sphericalAngularDamper::
sphericalAngularDamper
(
    const word& name,
    const dictionary& dict
)
:
    objectMotionRestraint(name, dict),
    coeff_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionRestraints::sphericalAngularDamper::
~sphericalAngularDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionRestraints::sphericalAngularDamper::restrain
(
    const movingObject& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintMoment = -coeff_*motion.omega();
    restraintForce = Zero;

    if (motion.report())
    {
        Info<< " moment " << restraintMoment
            << endl;
    }
}


bool Foam::objectMotionRestraints::sphericalAngularDamper::read
(
    const dictionary& dict
)
{
    objectMotionRestraint::read(dict);

    coeffDict_.lookup("coeff") >> coeff_;

    return true;
}


void Foam::objectMotionRestraints::sphericalAngularDamper::write
(
    Ostream& os
) const
{
    writeEntry(os, "coeff", coeff_);
}


// ************************************************************************* //
