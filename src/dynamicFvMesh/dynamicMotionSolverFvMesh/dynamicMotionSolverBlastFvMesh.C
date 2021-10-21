/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "dynamicMotionSolverBlastFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverBlastFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicBlastFvMesh,
        dynamicMotionSolverBlastFvMesh,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverBlastFvMesh::dynamicMotionSolverBlastFvMesh(const IOobject& io)
:
    dynamicBlastFvMesh(io),
    motionPtr_(motionSolver::New(*this, dynamicMeshDict())),
    velocityMotionCorrection_(*this, dynamicMeshDict())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverBlastFvMesh::~dynamicMotionSolverBlastFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::motionSolver& Foam::dynamicMotionSolverBlastFvMesh::motion() const
{
    return motionPtr_();
}


bool Foam::dynamicMotionSolverBlastFvMesh::refine(const bool)
{
    return false;
}


bool Foam::dynamicMotionSolverBlastFvMesh::update()
{
    fvMesh::movePoints(pointField(motionPtr_->newPoints()));
    velocityMotionCorrection_.update();

    return true;
}


bool Foam::dynamicMotionSolverBlastFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    motionPtr_->write();
    return fvMesh::writeObject(fmt, ver, cmp, write);
}


// ************************************************************************* //
