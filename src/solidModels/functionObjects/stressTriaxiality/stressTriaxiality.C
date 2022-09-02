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
    defineTypeNameAndDebug(stressTriaxiality, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        stressTriaxiality,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::stressTriaxiality::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volSymmTensorField& sigma =
            mesh_.lookupObject<volSymmTensorField>("sigma");

        // Calculate hydrostatic stress
        const volScalarField sigmaHyd(-tr(sigma)/3.0);

        // Calculate equivalent stress
        volScalarField sigmaEq(sqrt((3.0/2.0)*magSqr(dev(sigma))));

        // Limit sigmaEq to at least SMALL to avid division by zero
        sigmaEq = max(sigmaEq, dimensionedScalar("SMALL", dimPressure, SMALL));

        // Calculate stress triaxiality
        const volScalarField stressTriaxiality
        (
            "stressTriaxiality", -sigmaHyd/sigmaEq
        );

        stressTriaxiality.write();

        Info<< "Stress triaxiality : min = " << gMin(stressTriaxiality)
            << ", max = " << gMax(stressTriaxiality) << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stressTriaxiality::stressTriaxiality
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::stressTriaxiality::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::stressTriaxiality::execute()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::stressTriaxiality::read(const dictionary& dict)
{
    return true;
}


bool Foam::stressTriaxiality::write()
{
    return writeData();
}

// ************************************************************************* //
