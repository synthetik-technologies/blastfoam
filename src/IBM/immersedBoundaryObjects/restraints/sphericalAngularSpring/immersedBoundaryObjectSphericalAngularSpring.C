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

#include "immersedBoundaryObjectSphericalAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBoundaryObjectRestraints
{
    defineTypeNameAndDebug(sphericalAngularSpring, 0);

    addToRunTimeSelectionTable
    (
        immersedBoundaryObjectRestraint,
        sphericalAngularSpring,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectRestraints::sphericalAngularSpring::
sphericalAngularSpring
(
    const polyMesh& mesh,
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    immersedBoundaryObjectRestraint(mesh, name, sDoFRBMRDict),
    refQ_(),
    stiffness_(),
    damping_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectRestraints::sphericalAngularSpring::
~sphericalAngularSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::immersedBoundaryObjectRestraints::sphericalAngularSpring::restrain
(
    const immersedBoundaryObject& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintMoment = Zero;

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        vector axis = Zero;
        axis[cmpt] = 1;

        vector refDir = Zero;
        refDir[(cmpt + 1) % 3] = 1;

        vector newDir = motion.orientation() & refDir;

        axis = (refQ_ & axis);
        refDir = (refQ_ & refDir);
        newDir -= (axis & newDir)*axis;

        restraintMoment += -stiffness_*(refDir ^ newDir);
    }

    restraintMoment += -damping_*motion.omega();

    restraintForce = Zero;

    // Not needed to be altered as restraintForce is zero, but set to
    // centreOfRotation to be sure of no spurious moment
    restraintPosition = motion.centreOfRotation();

    if (motion.report())
    {
        Info<< " moment " << restraintMoment
            << endl;
    }
}


bool Foam::immersedBoundaryObjectRestraints::sphericalAngularSpring::read
(
    const dictionary& sDoFRBMRDict
)
{
    immersedBoundaryObjectRestraint::read(sDoFRBMRDict);

    refQ_ = sDoFRBMRCoeffs_.lookupOrDefault<tensor>("referenceOrientation", I);

    if (mag(mag(refQ_) - sqrt(3.0)) > 1e-9)
    {
        FatalErrorInFunction
            << "referenceOrientation " << refQ_ << " is not a rotation tensor. "
            << "mag(referenceOrientation) - sqrt(3) = "
            << mag(refQ_) - sqrt(3.0) << nl
            << exit(FatalError);
    }

    sDoFRBMRCoeffs_.lookup("stiffness") >> stiffness_;
    sDoFRBMRCoeffs_.lookup("damping") >> damping_;

    return true;
}


void Foam::immersedBoundaryObjectRestraints::sphericalAngularSpring::write
(
    Ostream& os
) const
{
    writeEntry(os, "referenceOrientation", refQ_);

    writeEntry(os, "stiffness", stiffness_);

    writeEntry(os, "damping", damping_);
}


// ************************************************************************* //
