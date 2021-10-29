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

#include "immersedBoundaryObjectLineConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBoundaryObjectConstraints
{
    defineTypeNameAndDebug(line, 0);

    addToRunTimeSelectionTable
    (
        immersedBoundaryObjectConstraint,
        line,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectConstraints::line::line
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const immersedBoundaryObject& motion
)
:
    immersedBoundaryObjectConstraint(name, sDoFRBMCDict, motion)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectConstraints::line::~line()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::immersedBoundaryObjectConstraints::line::setCentreOfRotation
(
    point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::immersedBoundaryObjectConstraints::line::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, direction_)));
}


void Foam::immersedBoundaryObjectConstraints::line::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::immersedBoundaryObjectConstraints::line::read
(
    const dictionary& sDoFRBMCDict
)
{
    immersedBoundaryObjectConstraint::read(sDoFRBMCDict);

    centreOfRotation_ = sDoFRBMCCoeffs_.lookupOrDefault
    (
        "centreOfRotation",
        motion_.initialCentreOfMass()
    );

    sDoFRBMCCoeffs_.lookup("direction") >> direction_;

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


void Foam::immersedBoundaryObjectConstraints::line::write
(
    Ostream& os
) const
{
    writeEntry(os, "centreOfRotation", centreOfRotation_);
    writeEntry(os, "direction", direction_);
}

// ************************************************************************* //
