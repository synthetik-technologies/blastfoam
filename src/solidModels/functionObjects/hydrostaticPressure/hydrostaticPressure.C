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

#include "hydrostaticPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hydrostaticPressure, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        hydrostaticPressure,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::hydrostaticPressure::writeData()
{
    if (time_.outputTime())
    {
        // Lookup the solid mesh
        const fvMesh* meshPtr = NULL;
        if (time_.foundObject<fvMesh>("solid"))
        {
            meshPtr = &(time_.lookupObject<fvMesh>("solid"));
        }
        else
        {
            meshPtr = &(time_.lookupObject<fvMesh>("region0"));
        }
        const fvMesh& mesh = *meshPtr;

        // Lookup stress tensor
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        // Calculate hydrostatic stress

        const volScalarField hydrostaticPressure
        (
            "hydrostaticPressure", -tr(sigma)/3.0
        );

        hydrostaticPressure.write();

        Info<< "Hydrostatic pressure: min = " << gMin(hydrostaticPressure)
            << ", max = " << gMax(hydrostaticPressure) << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrostaticPressure::hydrostaticPressure
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t)
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hydrostaticPressure::~hydrostaticPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::hydrostaticPressure::start()
{
    if (time_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::hydrostaticPressure::execute()
{
    if (time_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::hydrostaticPressure::read(const dictionary& dict)
{
    return true;
}


bool Foam::hydrostaticPressure::write()
{
    return writeData();
}

// ************************************************************************* //
