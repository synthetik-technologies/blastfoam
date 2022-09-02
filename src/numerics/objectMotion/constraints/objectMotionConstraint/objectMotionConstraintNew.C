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

#include "objectMotionConstraint.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::objectMotionConstraint>
Foam::objectMotionConstraint::New
(
    const word& name,
    const dictionary& dict,
    const movingObject& motion
)
{
    const word constraintType
    (
        dict.lookup("type")
    );

    motionConstructorTable::iterator cstrIter =
        motionConstructorTablePtr_->find(constraintType);

    if (cstrIter == motionConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown objectMotionConstraint type "
            << constraintType << nl << nl
            << "Valid objectMotionConstraints are : " << endl
            << motionConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<objectMotionConstraint>
    (
        cstrIter()(name, dict, motion)
    );
}


// ************************************************************************* //
