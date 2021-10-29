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

#include "immersedBoundaryObjectLinearSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBoundaryObjectRestraints
{
    defineTypeNameAndDebug(linearSpring, 0);

    addToRunTimeSelectionTable
    (
        immersedBoundaryObjectRestraint,
        linearSpring,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectRestraints::linearSpring::linearSpring
(
    const polyMesh& mesh,
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    immersedBoundaryObjectRestraint(mesh, name, sDoFRBMRDict),
    anchor_(),
    refAttachmentPt_(),
    stiffness_(),
    damping_(),
    restLength_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectRestraints::linearSpring::~linearSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::immersedBoundaryObjectRestraints::linearSpring::restrain
(
    const immersedBoundaryObject& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);

    vector r = restraintPosition - anchor_;

    scalar magR = mag(r);
    r /= (magR + vSmall);

    vector v = motion.velocity(restraintPosition);

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


bool Foam::immersedBoundaryObjectRestraints::linearSpring::read
(
    const dictionary& sDoFRBMRDict
)
{
    immersedBoundaryObjectRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.lookup("anchor") >> anchor_;
    sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;
    sDoFRBMRCoeffs_.lookup("stiffness") >> stiffness_;
    sDoFRBMRCoeffs_.lookup("damping") >> damping_;
    sDoFRBMRCoeffs_.lookup("restLength") >> restLength_;

    return true;
}


void Foam::immersedBoundaryObjectRestraints::linearSpring::write
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
