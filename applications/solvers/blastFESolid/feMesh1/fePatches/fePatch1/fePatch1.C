/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "fePatch1.H"
#include "addToRunTimeSelectionTable.H"
#include "feBoundaryMesh1.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fePatch1, 0);
    defineRunTimeSelectionTable(fePatch1, polyPatch);
    addToRunTimeSelectionTable(fePatch1, fePatch1, polyPatch);
}


Foam::autoPtr<Foam::fePatch1> Foam::fePatch1::New
(
    const polyPatch& patch,
    const feBoundaryMesh1& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fePatch1" << endl;
    }

    polyPatchConstructorTable::iterator cstrIter =
        polyPatchConstructorTablePtr_->find(patch.type());

    if (cstrIter == polyPatchConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fePatch1 type "
            << patch.type()
            << nl << nl
            << "Valid fePatch1 types are :" << endl
            << polyPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fePatch1>(cstrIter()(patch, bm));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fePatch1::fePatch1(const polyPatch& p, const feBoundaryMesh1& bm)
:
    polyPatch_(p),
    boundaryMesh_(bm),
    elementsPtr_(nullptr),
    localElementsPtr_(nullptr),
    localNodesPtr_(nullptr),
    shapesPtr_(nullptr),
    dshapesPtr_(nullptr),
    BsPtr_(nullptr),
    JsPtr_(nullptr),
    invJsPtr_(nullptr),
    WsPtr_(nullptr),
    nodeNormalsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fePatch1::~fePatch1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fePatch1::nNodes() const
{
    return meshNodes().size();
}


const Foam::elementList& Foam::fePatch1::elements() const
{
    if (!elementsPtr_)
    {
        calcElements();
    }
    return *elementsPtr_;
}


const Foam::List<Foam::labelList>& Foam::fePatch1::localElements() const
{
    if (!localElementsPtr_)
    {
        calcLocalElements();
    }
    return *localElementsPtr_;
}


const Foam::pointField& Foam::fePatch1::localNodes() const
{
    if (!localNodesPtr_)
    {
        calcLocalNodes();
    }
    return *localNodesPtr_;
}



const Foam::labelList& Foam::fePatch1::meshNodes() const
{
    return polyPatch_.meshPoints();
}


const Foam::List<Foam::List<Foam::scalarList>>& Foam::fePatch1::shapes() const
{
    if (!shapesPtr_)
    {
        calcShapes();
    }

    return *shapesPtr_;
}



const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::fePatch1::dshapes() const
{
    if (!dshapesPtr_)
    {
        calcDShapes();
    }

    return *dshapesPtr_;
}


const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::fePatch1::Bs() const
{
    if (!BsPtr_)
    {
        calcBs();
    }

    return *BsPtr_;
}


const Foam::List<Foam::List<Foam::tensor>>& Foam::fePatch1::Js() const
{
    if (!JsPtr_)
    {
        calcJs();
    }

    return *JsPtr_;
}


const Foam::List<Foam::List<Foam::tensor>>& Foam::fePatch1::invJs() const
{
    if (!invJsPtr_)
    {
        calcInvJs();
    }

    return *invJsPtr_;
}


const Foam::List<Foam::List<Foam::scalar>>& Foam::fePatch1::Ws() const
{
    if (!WsPtr_)
    {
        calcWs();
    }

    return *WsPtr_;
}


const Foam::List<Foam::List<Foam::vector>>&
Foam::fePatch1::nodeNormals() const
{
    if (!nodeNormalsPtr_)
    {
        calcNodeNormals();
    }
    return *nodeNormalsPtr_;
}


// ************************************************************************* //
