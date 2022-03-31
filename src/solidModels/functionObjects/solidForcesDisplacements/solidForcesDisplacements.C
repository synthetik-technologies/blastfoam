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

#include "solidForcesDisplacements.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidForcesDisplacements, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidForcesDisplacements,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidForcesDisplacements::writeData()
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

    // Lookup the displacement field
    const vectorField& D =
        mesh.lookupObject<volVectorField>
        (
            "D"
        ).boundaryField()[historyPatchID_];

    // Patch area vectors
    const vectorField& patchSf =
        mesh.Sf().boundaryField()[historyPatchID_];

    // Calculate the force as the intergal of the traction over the area
    vector force = vector::zero;

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
        const vectorField patchDeformSf(J*Finv.T() & patchSf);

        // It is assumed that sigma is the true (Cauchy) stress
        force = gSum(patchDeformSf & sigma);
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
    }
    else
    {
        // Linear geometry

        force = gSum(patchSf & sigma);
    }

    // Arithmetic average disp and force on patch
    vector avDisp = average(D);

    if (Pstream::master())
    {
        historyFilePtr_()
            << time_.time().value() << " "
                << avDisp.x() << " "
                << avDisp.y() << " "
                << avDisp.z() << " "
                << force.x() << " "
                << force.y() << " "
                << force.z();

        historyFilePtr_() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidForcesDisplacements::solidForcesDisplacements
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
    historyFilePtr_()
{
    Info << "Creating " << this->name() << " function object." << endl;

    word historyPatchName("notSpecified");
    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch") >> historyPatchName;
    }
    else
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "historyPatch not specified."
            << abort(FatalError);
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
        FatalErrorIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found."
            << abort(FatalError);
    }

    // Create history file if not already created
    if (historyFilePtr_.empty())
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
                    historyDir/"solidForcesDisplacements"
                  + historyPatchName + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "dispX" << " " << "dispY" << " "
                    << "dispZ" << " "
                    << "forceX" << " " << "forceY" << " "
                    << "forceZ";

                historyFilePtr_() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidForcesDisplacements::start()
{
    return writeData();
}


bool Foam::solidForcesDisplacements::execute()
{
    return writeData();
}


bool Foam::solidForcesDisplacements::read(const dictionary& dict)
{
    return true;
}


bool Foam::solidForcesDisplacements::write()
{
    return writeData();
}

// ************************************************************************* //
