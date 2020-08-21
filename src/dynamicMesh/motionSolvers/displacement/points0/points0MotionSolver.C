/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "points0MotionSolver.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(points0MotionSolver, 0);
}

// * * * * * * * * * * * * * * Private Data Members  * * * * * * * * * * * * //

Foam::IOobject Foam::points0MotionSolver::points0IO
(
    const polyMesh& mesh
) const
{
    const word instance0
    (
        mesh.time().findInstance
        (
            mesh.meshDir(),
            "points0",
            IOobject::READ_IF_PRESENT
        )
    );

    if (instance0 != mesh.time().constant())
    {
        // Points0 written to a time folder

        return IOobject
        (
            "points0",
            instance0,
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
    }
    else
    {
        // Check that points0 are actually in constant directory

        IOobject io
        (
            "points0",
            instance0,
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (io.typeHeaderOk<pointIOField>())
        {
            return io;
        }
        else
        {
            const word instance
            (
                mesh.time().findInstance
                (
                    mesh.meshDir(),
                    "points",
                    IOobject::READ_IF_PRESENT
                )
            );

            // Copy of original mesh points
            return IOobject
            (
                "points",
                instance,
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
        }
    }
}

const Foam::pointVectorField& Foam::points0MotionSolver::points0Init
(
    const polyMesh& mesh
) const
{
    pointIOField pio(points0IO(mesh));

    pointVectorField* p0Ptr
    (
        new pointVectorField
        (
            IOobject
            (
                "points0",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pointMesh::New(mesh),
            dimensionedVector(dimLength, Zero)
        )
    );
    p0Ptr->primitiveFieldRef() = pio;

    p0Ptr->store(p0Ptr);

    return mesh.lookupObject<pointVectorField>("points0");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::points0MotionSolver::points0MotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    points0_(points0Init(mesh))
{
    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file "
            <<  typeFilePath<pointIOField>
                (
                    IOobject
                    (
                        "points",
                        mesh.time().constant(),
                        polyMesh::meshSubDir,
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::points0MotionSolver::~points0MotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::points0MotionSolver::movePoints(const pointField&)
{}

void Foam::points0MotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // pointMesh already updates pointFields

    motionSolver::updateMesh(mpm);
}


// ************************************************************************* //
