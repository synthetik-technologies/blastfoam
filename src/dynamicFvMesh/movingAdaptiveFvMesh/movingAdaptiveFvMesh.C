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

#include "movingAdaptiveFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "displacementMotionSolver.H"
#include "componentDisplacementMotionSolver.H"
#include "MapGeometricFields.H"
#include "pointMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingAdaptiveFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicBlastFvMesh,
        movingAdaptiveFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingAdaptiveFvMesh::movingAdaptiveFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    adaptiveFvMesh(io),
    motionPtr_(motionSolver::New(*this, dynamicMeshDict())),
    velocityMotionCorrection_(*this, dynamicMeshDict())
{
    this->meshCutter().locationMapper().needMap();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingAdaptiveFvMesh::~movingAdaptiveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingAdaptiveFvMesh::updateMesh(const mapPolyMesh& mpm)
{
//     if (!isRefining_ && !isUnrefining_)
//     {
//         motionPtr_->updateMesh(mpm);
//         return;
//     }

    // Do not update while balancing this is handled in the
    // distribute function
    if (isBalancing_)
    {
//         return;
    }
    else if (isA<displacementMotionSolver>(motionPtr_()))
    {
        displacementMotionSolver& dispMS =
            dynamicCast<displacementMotionSolver>(motionPtr_());
        if (isRefining_)
        {
            meshCutter().locationMapper().interpolateMidPoints
            (
                dispMS.points0()
            );
        }
        else
        {
            pointMapper(dispMS.pointDisplacement().mesh(), mpm)
            (
                dispMS.points0(),
                dispMS.points0()
            );
        }

//         if (Pstream::parRun())
//         {
//             this->pushUntransformedData(dispMS.points0());
//         }
//         dispMS.pointDisplacement().primitiveFieldRef() =
//             this->points() - dispMS.points0();
    }
    else if (isA<componentDisplacementMotionSolver>(motionPtr_()))
    {
        componentDisplacementMotionSolver& dispMS =
            dynamicCast<componentDisplacementMotionSolver>
            (
                motionPtr_()
            );
        pointMapper(pointMesh::New(*this), mpm)(dispMS.points0());
//         if (Pstream::parRun())
//         {
//             this->pushUntransformedData(dispMS.points0());
//         }
    }

    adaptiveFvMesh::updateMesh(mpm);
}


void Foam::movingAdaptiveFvMesh::distribute
(
    const mapDistributePolyMesh& map
)
{
    adaptiveFvMesh::distribute(map);
    if (isA<displacementMotionSolver>(motionPtr_()))
    {
        displacementMotionSolver& dispMS =
            dynamicCast<displacementMotionSolver>(motionPtr_());
        map.distributePointData(dispMS.points0());
    }
    else if (isA<componentDisplacementMotionSolver>(motionPtr_()))
    {
        componentDisplacementMotionSolver& dispMS =
            dynamicCast<componentDisplacementMotionSolver>
            (
                motionPtr_()
            );
         map.distributePointData(dispMS.points0());
    }
}


const Foam::motionSolver& Foam::movingAdaptiveFvMesh::motion() const
{
    return motionPtr_();
}


bool Foam::movingAdaptiveFvMesh::refine()
{
    if (adaptiveFvMesh::refine())
    {
        meshCutter().locationMapper().clearOut();
        return true;
    }

    return false;
}


bool Foam::movingAdaptiveFvMesh::update()
{
    // Get the new points solving for displacement
    pointField pointsNew(motionPtr_->newPoints());

    //- Sync points across boundaries
//     if (Pstream::parRun())
//     {
//        this->pushUntransformedData(pointsNew);
//     }

    //- Move mesh
    fvMesh::movePoints(pointsNew);
    velocityMotionCorrection_.update();


    return moving();
}


bool Foam::movingAdaptiveFvMesh::writeObject
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
        pointIOField points0
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
        );
        pointVectorField ppoints0
        (
            IOobject
            (
                "points0",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(*this),
            dimensionedVector(dimLength, Zero)
        );
        ppoints0.primitiveFieldRef() = points0;
        ppoints0.write();
        points0.write();


    }
    return adaptiveFvMesh::writeObject(fmt, ver, cmp, write);
}


// ************************************************************************* //
