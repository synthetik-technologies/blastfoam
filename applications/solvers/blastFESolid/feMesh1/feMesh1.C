/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "feMesh1.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(feMesh1, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::feMesh1::clearRefGeom()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "clearRefGeom" << endl;
    }

    deleteDemandDrivenData(shapesPtr_);
    deleteDemandDrivenData(dshapesPtr_);

    clearPhysGeom();
}


void Foam::feMesh1::clearPhysGeom()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "clearPhysGeom" << endl;
    }

    deleteDemandDrivenData(BsPtr_);
    deleteDemandDrivenData(JsPtr_);
    deleteDemandDrivenData(invJsPtr_);
    deleteDemandDrivenData(WsPtr_);
    deleteDemandDrivenData(WPtr_);
}


void Foam::feMesh1::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "isMeshUpdate: " << isMeshUpdate << endl;
    }

    deleteDemandDrivenData(elementsPtr_);
    deleteDemandDrivenData(ipLabelsPtr_);
    nIp_ = -1;
    clearRefGeom();
}


void Foam::feMesh1::clearOut()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::feMesh1::feMesh1(const polyMesh& pMesh, const label intOrder)
:
    MeshObject<polyMesh, PatchMeshObject, feMesh1>(pMesh),

    pointMesh_(pointMesh::New(pMesh)),
    intOrder_(intOrder),
    shellIntOrder_(-1),
    shell_(false),
    boundary_(*this, pMesh.boundaryMesh()),

    elementsPtr_(nullptr),
    shapesPtr_(nullptr),
    dshapesPtr_(nullptr),
    BsPtr_(nullptr),
    JsPtr_(nullptr),
    invJsPtr_(nullptr),
    WsPtr_(nullptr),
    WPtr_(nullptr),
    nIp_(-1),
    ipLabelsPtr_(nullptr)
{}


Foam::feMesh1::feMesh1
(
    const polyMesh& pMesh,
    const label orderRS,
    const label orderT
)
:
    MeshObject<polyMesh, PatchMeshObject, feMesh1>(pMesh),

    pointMesh_(pointMesh::New(pMesh)),
    intOrder_(orderRS),
    shellIntOrder_(orderT),
    shell_(true),
    boundary_(*this, pMesh.boundaryMesh()),

    elementsPtr_(nullptr),
    shapesPtr_(nullptr),
    dshapesPtr_(nullptr),
    BsPtr_(nullptr),
    JsPtr_(nullptr),
    invJsPtr_(nullptr),
    WsPtr_(nullptr),
    WPtr_(nullptr),
    nIp_(-1),
    ipLabelsPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::feMesh1::~feMesh1()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::feMesh1::reset(const bool validBoundary)
{
    const polyMesh& pm = mesh();
    if (debug)
    {
        Pout<< "feMesh1::reset(const bool validBoundary): "
            << "Resetting from polyMesh " << pm.name() << endl;
    }

    clearAddressing();

    boundary_.reset(pm.boundaryMesh());
    if (validBoundary)
    {
        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();
    }
}


bool Foam::feMesh1::movePoints()
{
    if (debug)
    {
        Pout<< "feMesh1::movePoints(const pointField&): "
            << "Moving points." << endl;
    }

    clearPhysGeom();

    boundary_.movePoints(mesh().points());

    return true;
}


void Foam::feMesh1::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "feMesh1::updateMesh(const mapPolyMesh&): "
            << "Updating for topology changes." << endl;
        Pout<< endl;
    }

    clearAddressing();

    boundary_.updateMesh();
}


void Foam::feMesh1::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    if (debug)
    {
        Pout<< "feMesh1::reorderPatches( const labelUList&, const bool): "
            << "Updating for reordered patches." << endl;
        Pout<< endl;
    }

    boundary_.shuffle(newToOld, validBoundary);
}


void Foam::feMesh1::addPatch(const label patchi)
{
    if (debug)
    {
        Pout<< "feMesh1::addPatch(const label): "
            << "Adding patch at " << patchi << endl;
        Pout<< endl;
    }

    boundary_.set
    (
        patchi,
        fePatch1::New
        (
            mesh().boundaryMesh()[patchi], boundary_
        ).ptr()
     );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::feMesh1::operator!=(const feMesh1& bm) const
{
    return &bm != this;
}


bool Foam::feMesh1::operator==(const feMesh1& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
