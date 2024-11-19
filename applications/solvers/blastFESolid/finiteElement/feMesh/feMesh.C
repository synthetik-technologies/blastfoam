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

#include "feMesh.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"
#include "mapClouds.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(feMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::feMesh::mapFields(const mapPolyMesh& mpm)
{
//     if (debug)
//     {
//         Pout<< "void feMesh::mapFields(const mapPolyMesh&): "
//             << "Mapping all registered feFields."
//             << endl;
//     }
//     // Create a mapper
//     const feMeshMapper m(*this, mpm);
//
//     MapGeometricFields<scalar, fePatchField, feMeshMapper, feMesh>(m);
//     MapGeometricFields<vector, fePatchField, feMeshMapper, feMesh>(m);
//     MapGeometricFields
//     <
//         sphericalTensor,
//         fePatchField,
//         feMeshMapper,
//         feMesh
//     >(m);
//     MapGeometricFields<symmTensor, fePatchField, feMeshMapper, feMesh>
//     (m);
//     MapGeometricFields<tensor, fePatchField, feMeshMapper, feMesh>(m);
}


void Foam::feMesh::clearGeom()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "clearGeom" << endl;
    }

//     meshObject::clearUpto
//     <
//         feMesh,
//         GeometricMeshObject,
//         MoveableMeshObject
//     >(*this);
//
//
//     deleteDemandDrivenData(shapesPtr_);
//     deleteDemandDrivenData(dshapesPtr_);
//     deleteDemandDrivenData(BsPtr_);
//     deleteDemandDrivenData(JsPtr_);
//     deleteDemandDrivenData(invJsPtr_);
//     deleteDemandDrivenData(WsPtr_);
}


void Foam::feMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "isMeshUpdate: " << isMeshUpdate << endl;
    }

//     if (isMeshUpdate)
//     {
//         // Part of a mesh update. Keep meshObjects that have an updateMesh
//         // callback
//         meshObject::clearUpto
//         <
//             feMesh,
//             TopologicalMeshObject,
//             UpdateableMeshObject
//         >
//         (
//             *this
//         );
//     }
//     else
//     {
//         meshObject::clear<feMesh, TopologicalMeshObject>(*this);
//     }
    deleteDemandDrivenData(elementsPtr_);
    deleteDemandDrivenData(edgeNodesPtr_);
    deleteDemandDrivenData(faceNodesPtr_);
    deleteDemandDrivenData(nodesPtr_);
    deleteDemandDrivenData(ipLabelsPtr_);
    nIp_ = -1;
}


void Foam::feMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::feMesh::feMesh(const label order, const polyMesh& pMesh)
:
    objectRegistry
    (
        IOobject
        (
            IOobject::groupName("feMesh", Foam::name(order)),
            pMesh.facesInstance(),
            pMesh
        )
    ),
    MeshObject<polyMesh, Foam::PatchMeshObject, feMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),

    order_(order),
    boundary_(*this, pMesh.boundaryMesh()),

    elementsPtr_(nullptr),
    nodesPtr_(nullptr),
    edgeNodesPtr_(nullptr),
    faceNodesPtr_(nullptr),
    shapesPtr_(nullptr),
    dshapesPtr_(nullptr),
    BsPtr_(nullptr),
    JsPtr_(nullptr),
    invJsPtr_(nullptr),
    WsPtr_(nullptr),
    nIp_(-1),
    ipLabelsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::feMesh::~feMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::feMesh::reset(const bool validBoundary)
{
    const polyMesh& pm = mesh();
    if (debug)
    {
        Pout<< "feMesh::reset(const bool validBoundary): "
            << "Resetting from polyMesh " << pm.name() << endl;
    }

    boundary_.reset(pm.boundaryMesh());
    if (validBoundary)
    {
        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();
    }
}


bool Foam::feMesh::movePoints()
{
    if (debug)
    {
        Pout<< "feMesh::movePoints(const pointField&): "
            << "Moving points." << endl;
    }

    boundary_.movePoints(GeoMesh<polyMesh>::mesh_.points());

    return true;
}


void Foam::feMesh::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "feMesh::updateMesh(const mapPolyMesh&): "
            << "Updating for topology changes." << endl;
        Pout<< endl;
    }
    boundary_.updateMesh();

    // Map all registered point fields
    mapFields(mpm);
}


void Foam::feMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    if (debug)
    {
        Pout<< "feMesh::reorderPatches( const labelUList&, const bool): "
            << "Updating for reordered patches." << endl;
        Pout<< endl;
    }

    boundary_.shuffle(newToOld, validBoundary);

//     objectRegistry& db = const_cast<objectRegistry&>(thisDb());
//     ReorderPatchFields<feScalarField>(db, newToOld);
//     ReorderPatchFields<feVectorField>(db, newToOld);
//     ReorderPatchFields<feSphericalTensorField>(db, newToOld);
//     ReorderPatchFields<feSymmTensorField>(db, newToOld);
//     ReorderPatchFields<feTensorField>(db, newToOld);
}


void Foam::feMesh::addPatch(const label patchi)
{
    if (debug)
    {
        Pout<< "feMesh::addPatch(const label): "
            << "Adding patch at " << patchi << endl;
        Pout<< endl;
    }

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    if (pbm.size() != boundary_.size())
    {
        FatalErrorInFunction << "Problem :"
            << " feBoundaryMesh size :" << boundary_.size()
            << " polyBoundaryMesh size :" << pbm.size()
            << exit(FatalError);
    }

    boundary_.set(patchi, fePatch::New(pbm[patchi], boundary_).ptr());

//     objectRegistry& db = const_cast<objectRegistry&>(thisDb());
//     const dictionary d;
//     const word patchFieldType("calculated");

//     AddPatchFields<feScalarField>(db, patchi, d, patchFieldType, Zero);
//     AddPatchFields<feVectorField>(db, patchi, d, patchFieldType, Zero);
//     AddPatchFields<feSphericalTensorField>
//     (
//         db,
//         patchi,
//         d,
//         patchFieldType,
//         Zero
//     );
//     AddPatchFields<feSymmTensorField>(db, patchi, d, patchFieldType, Zero);
//     AddPatchFields<feTensorField>(db, patchi, d, patchFieldType, Zero);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::feMesh::operator!=(const feMesh& bm) const
{
    return &bm != this;
}


bool Foam::feMesh::operator==(const feMesh& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
