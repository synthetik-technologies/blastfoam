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

#include "solidKineticEnergy.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidKineticEnergy, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidKineticEnergy,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidKineticEnergy::writeData()
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

    // Skip if the velocity field is not found
    if (!mesh.foundObject<volVectorField>("U"))
    {
        return false;
    }

    // Lookup the velocity field
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    // Look rho from the solid model
    const volScalarField& rho = lookupSolidModel(mesh).rho();

    // Calculate the kinetic energy per unit volume field
    const volScalarField kinEnergyPerVol = 0.5*rho*(U & U);

    // Calculate the total kinetic energy
    const scalar kinEnergy =
        gSum(kinEnergyPerVol.internalField()*mesh.V().field());

    // Write to file
    if (Pstream::master())
    {
        historyFilePtr_()
            << time_.time().value() << " " << kinEnergy << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidKineticEnergy::solidKineticEnergy
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object." << endl;

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            const word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir/"solidKineticEnergy.dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " " << "kineticEnergy" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidKineticEnergy::start()
{
    return writeData();
}


bool Foam::solidKineticEnergy::execute()
{
    return writeData();
}


bool Foam::solidKineticEnergy::read(const dictionary& dict)
{
    return true;
}


bool Foam::solidKineticEnergy::write()
{
    return writeData();
}

// ************************************************************************* //
