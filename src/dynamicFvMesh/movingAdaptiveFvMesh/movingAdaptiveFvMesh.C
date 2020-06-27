/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020   Jeff Heylmun:  Set old cell volumes after refinement is done.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingAdaptiveFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        movingAdaptiveFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingAdaptiveFvMesh::movingAdaptiveFvMesh(const IOobject& io)
:
    adaptiveFvMesh(io),
    motionPtr_(motionSolver::New(*this, dynamicMeshDict())),
    velocityMotionCorrection_(*this, dynamicMeshDict())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingAdaptiveFvMesh::~movingAdaptiveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingAdaptiveFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    fvMesh::updateMesh(mpm);
    motionPtr_->updateMesh(mpm);
}

const Foam::motionSolver& Foam::movingAdaptiveFvMesh::motion() const
{
    return motionPtr_();
}


bool Foam::movingAdaptiveFvMesh::update()
{
    // Refine mesh
    adaptiveFvMesh::update();
    volScalarField::Internal Vold(this->V());

    //- Move mesh
    fvMesh::movePoints(motionPtr_->newPoints());
    velocityMotionCorrection_.update();

    this->setV0() = Vold;

    return true;
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
    return adaptiveFvMesh::writeObject(fmt, ver, cmp, write);
}


// ************************************************************************* //
