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

#include "solidPotentialEnergy.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPotentialEnergy, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPotentialEnergy,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPotentialEnergy::writeData()
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

    // Lookup the solid model
    const solidModel& solMod = lookupSolidModel(mesh);

    // Only currently implemented for linear geometry solid models
    if (solMod.nonLinGeom() != nonLinearGeometry::LINEAR_GEOMETRY)
    {
        WarningIn(type())
            << "Only implemented for linear geometry solid models" << endl;
        return false;
    }

    // Lookup the gravity field from the solid model
    const uniformDimensionedVectorField& g = solMod.g();

    // Look rho from the solid model
    const volScalarField& rho = solMod.rho();

    // Calculate the unit direction of gravity
    const scalar magG = mag(g.value());
    if (magG < SMALL)
    {
        WarningIn(type())
            << "Gravity is zero: skipping" << endl;
        return false;
    }
    const vector gDir = g.value()/magG;

    // Lookup the displacement field
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    // Calculate the height from the reference (zero potential energy) plane for
    // all cells
    // Potential energy is positive for any cells that are in the negative
    // gravity direction from the refPoint
    const scalarField h =
        (-gDir & (mesh.C().internalField() + U.internalField() - refPoint_));

    // Calculate the potential energy per unit volume field
    const scalarField potentialEnergyPerVol = rho.internalField()*magG*h;

    // Calculate the total potential energy
    const scalar potentialEnergy = gSum(potentialEnergyPerVol*mesh.V().field());

    // Write to file
    if (Pstream::master())
    {
        historyFilePtr_()
            << time_.time().value() << " " << potentialEnergy << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPotentialEnergy::solidPotentialEnergy
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    refPoint_(dict.lookup("referencePoint")),
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
                    historyDir/"solidPotentialEnergy.dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " " << "potentialEnergy" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPotentialEnergy::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::solidPotentialEnergy::execute(const bool forceWrite)
#else
bool Foam::solidPotentialEnergy::execute()
#endif
{
    return writeData();
}


bool Foam::solidPotentialEnergy::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::solidPotentialEnergy::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
