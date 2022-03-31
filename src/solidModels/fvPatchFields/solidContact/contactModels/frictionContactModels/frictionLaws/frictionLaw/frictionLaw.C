/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

InClass
    frictionLaw

\*---------------------------------------------------------------------------*/

#include "frictionLaw.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(frictionLaw, 0);
defineRunTimeSelectionTable(frictionLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
frictionLaw::frictionLaw
(
    const word& name,
    const frictionContactModel& frictionContactModel,
    const dictionary& dict
)
:
    frictionModel_(frictionContactModel),
    frictionLawDict_(dict)
{}


// Construct as a copy
frictionLaw::frictionLaw(const frictionLaw& fricLaw)
:
    frictionModel_(fricLaw.frictionModel_),
    frictionLawDict_(fricLaw.frictionLawDict_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar frictionLaw::slipTraction(const scalar pressure)
{
    notImplemented("scalar slipTraction(const scalar pressure)");

    // Keep compiler happy
    return 0.0;
}

scalar frictionLaw::slipTraction(const scalar pressure, const vector& slipDir)
{
    // If this function is not overwritten, then we default to the
    // more simple law
    return slipTraction(pressure);
}

scalar frictionLaw::slipTraction
(
        const scalar contactPressure,
        const vector& faceSlip,
        const vector& slaveFaceVelocity,
        const vector& masterFaceVelocity
)
{
    // If this function is not overwritten, then we default to the
    // more simple law
    return slipTraction(contactPressure);
}

scalar frictionLaw::slipTraction
(
    const scalar contactPressure,         // Contact pressure
    const vector& faceSlip,               // Slip vector
    const vector& slaveFaceVelocity,      // Velocity of slave face
    const vector& masterFaceVelocity,     // Velocity of master face
    const label slavePatchIndex,          // Slave patch index
    const label faceIndex                 // Local slave face ID
)
{
    // If this function is not overwritten, then we default to the
    // more simple law
    return slipTraction(contactPressure, faceSlip);
}

} // End namespace Foam

// ************************************************************************* //
