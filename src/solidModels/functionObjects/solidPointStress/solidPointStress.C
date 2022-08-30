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

#include "solidPointStress.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPointStress, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPointStress,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPointStress::writeData()
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

    if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        // Read the stress field
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        // Create a point mesh
        pointMesh pMesh(mesh);

        // Create a point stress field
        pointSymmTensorField pointSigma
        (
            IOobject
            (
                "pointSigma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("zero", sigma.dimensions(), symmTensor::zero)
        );

        mesh.lookupObject<volPointInterpolation>
        (
            "volPointInterpolation"
        ).interpolate(sigma, pointSigma);

        symmTensor pointSigmaValue = symmTensor::zero;
        if (pointID_ > -1)
        {
            pointSigmaValue = pointSigma[pointID_];
        }
        reduce(pointSigmaValue, sumOp<symmTensor>());

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                << " " << pointSigmaValue.xx()
                << " " << pointSigmaValue.xy()
                << " " << pointSigmaValue.xz()
                << " " << pointSigmaValue.yy()
                << " " << pointSigmaValue.yz()
                << " " << pointSigmaValue.zz()
                << endl;
        }
    }
    else
    {
        InfoIn(this->name() + " function object constructor")
            << "volSymmTensorField sigma not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPointStress::solidPointStress
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    pointID_(-1),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    // Lookup the point
    const vector point(dict.lookup("point"));

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

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // Find the closest point
        scalar minDist = GREAT;

        forAll(mesh.points(), pI)
        {
            scalar dist = mag(mesh.points()[pI] - point);

            if (dist < minDist)
            {
                minDist = dist;
                pointID_ = pI;
            }
        }

        // Find global closest point
        const scalar globalMinDist = returnReduce(minDist, minOp<scalar>());
        int procNo = -1;
        if (mag(globalMinDist - minDist) < SMALL)
        {
            procNo = Pstream::myProcNo();
        }
        else
        {
            pointID_ = -1;
        }

        // More than one processor can have the point so we will take the proc
        // with the lowest processor number
        const int globalMinProc = returnReduce(procNo, minOp<int>());
        if (mag(globalMinProc - procNo) > SMALL)
        {
            pointID_ = -1;
        }

        if (pointID_ > -1)
        {
            Pout<< this->name()
                << ": distance from specified point is " << minDist
                << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            const word startTimeName =
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
                new OFstream(historyDir/"solidPointStress_" + name + ".dat")
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time XX XY XZ YY YZ ZZ" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPointStress::start()
{
    return false;
}


bool Foam::solidPointStress::execute()
{
    return writeData();
}


bool Foam::solidPointStress::read(const dictionary& dict)
{
    return true;
}


bool Foam::solidPointStress::write()
{
    return writeData();
}

// ************************************************************************* //
