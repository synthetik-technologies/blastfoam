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

#include "solidForces.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidForces, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidForces,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidForces::writeData()
{
    if (patchFound_)
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

        // Patch area vectors
        const vectorField& patchSf =
            mesh.Sf().boundaryField()[historyPatchID_];

        // Patch unit area vectors
        const vectorField patchNf(mesh.boundary()[historyPatchID_].nf());

        // Calculate the force as the intergal of the traction over the area
        vector force = vector::zero;
        scalar normalForce = 0.0;

        // Lookup the stress field
        const symmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[historyPatchID_];

        // Check if it is a linear or nonlinear geometry case
        if
        (
            !mesh.foundObject<volVectorField>("DD")
         && mesh.foundObject<volTensorField>("F")
        )
        {
            // Total Lagrangian

            // Lookup the inverse of the deformation gradient
            const tensorField& Finv =
                mesh.lookupObject<volTensorField>
                (
                    "Finv"
                ).boundaryField()[historyPatchID_];

            // Lookup the Jacobian
            const scalarField& J =
                mesh.lookupObject<volScalarField>
                (
                    "J"
                ).boundaryField()[historyPatchID_];

            // Calculate area vectors in the deformed configuration
            const vectorField patchDeformSf((J*Finv.T() & patchSf));

            // Calculate unit area vectors in the deformed configuration
            const vectorField patchDeformNf(patchDeformSf/mag(patchDeformSf));

            // It is assumed that sigma is the true (Cauchy) stress
            force = gSum(patchDeformSf & sigma);

            normalForce = gSum(patchDeformNf & (patchDeformSf & sigma));
        }
        else if
        (
            mesh.foundObject<volVectorField>("DD")
         && mesh.foundObject<volTensorField>("F")
        )
        {
            // Updated Lagrangian

            // The mesh will be in the deformed configuration at the end of
            // the time-step

            // It is assumed that sigma is the true (Cauchy) stress
            force = gSum(patchSf & sigma);

            normalForce = gSum(patchNf & (patchSf & sigma));
        }
        else
        {
            // Linear geometry

            force = gSum(patchSf & sigma);

            normalForce = gSum(patchNf & (patchSf & sigma));
        }

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                    << " " << force.x()
                    << " " << force.y()
                    << " " << force.z()
                    << " " << normalForce
                    << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidForces::solidForces
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyPatchID_(-1),
    patchFound_(false),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    word historyPatchName("notSpecified");

    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch")
            >> historyPatchName;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "solidForces: historyPatch not specified" << endl;
    }

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

    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);

    if (historyPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found"
            << endl;
    }
    else
    {
        patchFound_ = true;
    }

    // Create history file if not already created
    if (historyFilePtr_.empty() && patchFound_)
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());


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
                        historyDir/"solidForces" + historyPatchName + ".dat"
                    )
                );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time forceX forceY forceZ normalForce" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidForces::start()
{
    return writeData();
}


bool Foam::solidForces::execute()
{
    return writeData();
}


bool Foam::solidForces::read(const dictionary& dict)
{
    return true;
}


bool Foam::solidForces::write()
{
    return writeData();
}

// ************************************************************************* //
