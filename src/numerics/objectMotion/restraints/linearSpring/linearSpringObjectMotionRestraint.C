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

#include "linearSpringObjectMotionRestraint.H"
#include "addToRunTimeSelectionTable.H"
#include "movingObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace objectMotionRestraints
{
    defineTypeNameAndDebug(linearSpring, 0);

    addToRunTimeSelectionTable
    (
        objectMotionRestraint,
        linearSpring,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectMotionRestraints::linearSpring::linearSpring
(
    const word& name,
    const dictionary& dict
)
:
    objectMotionRestraint( name, dict),
    anchor_(),
    refAttachmentPt_(),
    stiffness_(),
    damping_(),
    restLength_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectMotionRestraints::linearSpring::~linearSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionRestraints::linearSpring::restrain
(
    const movingObject& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);

    vector r = restraintPosition - anchor_;

    scalar magR = mag(r);
    r /= (magR + vSmall);

    vector v = motion.v(restraintPosition);

    restraintForce = -stiffness_*(magR - restLength_)*r - damping_*(r & v)*r;

    restraintMoment = Zero;

    if (motion.report())
    {
        Info<< " attachmentPt - anchor " << r*magR
            << " spring length " << magR
            << " force " << restraintForce
            << endl;
    }
}


bool Foam::objectMotionRestraints::linearSpring::read
(
    const dictionary& dict
)
{
    objectMotionRestraint::read(dict);

    coeffDict_.lookup("anchor") >> anchor_;
    coeffDict_.lookup("refAttachmentPt") >> refAttachmentPt_;
    coeffDict_.lookup("stiffness") >> stiffness_;
    coeffDict_.lookup("damping") >> damping_;
    coeffDict_.lookup("restLength") >> restLength_;

    return true;
}


void Foam::objectMotionRestraints::linearSpring::write
(
    Ostream& os
) const
{
    writeEntry(os, "anchor", anchor_);

    writeEntry(os, "refAttachmentPt", refAttachmentPt_);

    writeEntry(os, "stiffness", stiffness_);

    writeEntry(os, "damping", damping_);

    writeEntry(os, "restLength", restLength_);
}

// ************************************************************************* //
