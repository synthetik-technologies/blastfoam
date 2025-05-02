/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "stressTriaxiality.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stressTriaxiality, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        stressTriaxiality,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stressTriaxiality::stressTriaxiality
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stressTriaxiality::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::stressTriaxiality::execute()
{
    // Lookup stress tensor
    const volSymmTensorField& sigma =
        mesh_.lookupObject<volSymmTensorField>("sigma");

    // Calculate hydrostatic stress
    const volScalarField sigmaHyd(-tr(sigma)/3.0);

    // Calculate equivalent stress
    volScalarField sigmaEq(sqrt((3.0/2.0)*magSqr(dev(sigma))));

    // Limit sigmaEq to at least small to avid division by zero
    sigmaEq.max(dimensionedScalar(dimPressure, small));

    // Calculate stress triaxiality
    return store("stressTriaxiality", -sigmaHyd/sigmaEq);
}


bool Foam::functionObjects::stressTriaxiality::write()
{
    return writeObject("stressTriaxiality");
}

// ************************************************************************* //
