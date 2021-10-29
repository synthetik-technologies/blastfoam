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

#include "immersedBoundaryObjectConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryObjectConstraint, 0);
    defineRunTimeSelectionTable(immersedBoundaryObjectConstraint, motion);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectConstraint::immersedBoundaryObjectConstraint
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const immersedBoundaryObject& motion
)
:
    name_(name),
    sDoFRBMCCoeffs_(sDoFRBMCDict),
    motion_(motion)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectConstraint::~immersedBoundaryObjectConstraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::immersedBoundaryObjectConstraint::read
(
    const dictionary& sDoFRBMCDict
)
{
    sDoFRBMCCoeffs_ = sDoFRBMCDict;

    return true;
}


void Foam::immersedBoundaryObjectConstraint::write(Ostream& os) const
{}


// ************************************************************************* //
