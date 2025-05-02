/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "fePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "feBoundaryMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fePatch, 0);
    defineRunTimeSelectionTable(fePatch, polyPatch);
    addToRunTimeSelectionTable(fePatch, fePatch, polyPatch);
}


Foam::autoPtr<Foam::fePatch> Foam::fePatch::New
(
    const polyPatch& patch,
    const feBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fePatch" << endl;
    }

    polyPatchConstructorTable::iterator cstrIter =
        polyPatchConstructorTablePtr_->find(patch.type());

    if (cstrIter == polyPatchConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fePatch type "
            << patch.type()
            << nl << nl
            << "Valid fePatch types are :" << endl
            << polyPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fePatch>(cstrIter()(patch, bm));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fePatch::fePatch(const polyPatch& p, const feBoundaryMesh& bm)
:
    polyPatch_(p),
    boundaryMesh_(bm),
    elementsPtr_(nullptr),
    meshNodesPtr_(nullptr),
    shapesPtr_(nullptr),
    dshapesPtr_(nullptr),
    BsPtr_(nullptr),
    JsPtr_(nullptr),
    invJsPtr_(nullptr),
    WsPtr_(nullptr),
    nodeNormalsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fePatch::~fePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fePatch::nNodes() const
{
    return meshNodes().size();
}


const Foam::elementList& Foam::fePatch::elements() const
{
    if (!elementsPtr_)
    {
        calcElements();
    }
    return *elementsPtr_;
}


const Foam::labelList& Foam::fePatch::meshNodes() const
{
    if (!meshNodesPtr_)
    {
        calcMeshData();
    }
    return *meshNodesPtr_;
}


const Foam::List<Foam::List<Foam::scalarList>>& Foam::fePatch::shapes() const
{
    if (!shapesPtr_)
    {
        calcShapes();
    }

    return *shapesPtr_;
}



const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::fePatch::dshapes() const
{
    if (!dshapesPtr_)
    {
        calcDShapes();
    }

    return *dshapesPtr_;
}


const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::fePatch::Bs() const
{
    if (!BsPtr_)
    {
        calcBs();
    }

    return *BsPtr_;
}


const Foam::List<Foam::List<Foam::tensor>>& Foam::fePatch::Js() const
{
    if (!JsPtr_)
    {
        calcJs();
    }

    return *JsPtr_;
}


const Foam::List<Foam::List<Foam::tensor>>& Foam::fePatch::invJs() const
{
    if (!invJsPtr_)
    {
        calcInvJs();
    }

    return *invJsPtr_;
}


const Foam::List<Foam::List<Foam::scalar>>& Foam::fePatch::Ws() const
{
    if (!WsPtr_)
    {
        calcWs();
    }

    return *WsPtr_;
}


const Foam::List<Foam::List<Foam::vector>>&
Foam::fePatch::nodeNormals() const
{
    if (!nodeNormalsPtr_)
    {
        calcNodeNormals();
    }
    return *nodeNormalsPtr_;
}


// ************************************************************************* //
