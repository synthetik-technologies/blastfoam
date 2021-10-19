/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020   21-05-2020 Synthetik Applied Technologies: |  Set old cell volumes
                                after refinement is done.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "movingAdaptiveBlastFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "displacementMotionSolver.H"
#include "componentDisplacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingAdaptiveBlastFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicBlastFvMesh,
        movingAdaptiveBlastFvMesh,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingAdaptiveBlastFvMesh::movingAdaptiveBlastFvMesh(const IOobject& io)
:
    adaptiveBlastFvMesh(io),
    motionPtr_(motionSolver::New(*this, dynamicMeshDict())),
    velocityMotionCorrection_(*this, dynamicMeshDict())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingAdaptiveBlastFvMesh::~movingAdaptiveBlastFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingAdaptiveBlastFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    fvMesh::updateMesh(mpm);

    // Do not update while balancing since the point0 class does not
    // support this
    if (!isBalancing_)
    {
        motionPtr_->updateMesh(mpm);
    }
}

const Foam::motionSolver& Foam::movingAdaptiveBlastFvMesh::motion() const
{
    return motionPtr_();
}


bool Foam::movingAdaptiveBlastFvMesh::refine(const bool correctError)
{
    if (adaptiveBlastFvMesh::refine(correctError))
    {
        if (isA<displacementMotionSolver>(motionPtr_()))
        {
            displacementMotionSolver& dispMS =
                dynamicCast<displacementMotionSolver>(motionPtr_());
            pointField newPoints0
            (
                points() - dispMS.pointDisplacement().primitiveField()
            );
            this->pushUntransformedData(newPoints0);
            dispMS.points0() = newPoints0;
        }
        else if (isA<componentDisplacementMotionSolver>(motionPtr_()))
        {
            NotImplemented;
//             componentDisplacementMotionSolver& dispMS =
//                 dynamicCast<componentDisplacementMotionSolver>
//                 (
//                     motionPtr_()
//                 );
//             dispMS.points0() =
//                 points().component(dispMS.cmpt())
//               - dispMS.pointDisplacement().primitiveField();
        }
        return true;
    }

    return false;
}


bool Foam::movingAdaptiveBlastFvMesh::update()
{
    // Get the new points solving for displacement
    pointField pointsNew(motionPtr_->newPoints());

    //- Sync points across boundaries
    if (Pstream::parRun())
    {
       this->pushUntransformedData(pointsNew);
    }

    //- Move mesh
    fvMesh::movePoints(pointsNew);
    velocityMotionCorrection_.update();

    return moving();
}


bool Foam::movingAdaptiveBlastFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    motionPtr_->write();
    if (isA<displacementMotionSolver>(motionPtr_()))
    {
        const displacementMotionSolver& dispMS =
            dynamicCast<const displacementMotionSolver>(motionPtr_());
        pointIOField
        (
            IOobject
            (
                "points0",
                this->time().timeName(),
                polyMesh::meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dispMS.points0()
        ).write();
    }
    return adaptiveBlastFvMesh::writeObject(fmt, ver, cmp, write);
}


// ************************************************************************* //
