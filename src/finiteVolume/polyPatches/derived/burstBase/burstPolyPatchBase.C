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

#include "burstPolyPatchBase.H"
#include "globalMeshData.H"
#include "noBurstModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstPolyPatchBase, 0);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::burstPolyPatchBase::makePointIntact() const
{
    const polyBoundaryMesh& pbm = patch_.boundaryMesh();
    const polyMesh& mesh = pbm.mesh();

    // Set the point intact field
    pointIntact_.set(new scalarField(patch_.nPoints(), 0));
    scalarField& pI = pointIntact_();
    const scalarField& I = intact_();

    const labelListList& pointFaces = patch_.pointFaces();
    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];
        forAll(pFaces, fi)
        {
            const label facei = pFaces[fi];
            if (I[facei] > 0.5)
            {
                pI[pointi] = 1.0;
                break;
            }
        }
    }

    if (Pstream::parRun())
    {
        scalarField allPI(mesh.nPoints(), 0);
        UIndirectList<scalar>(allPI, patch_.meshPoints()) = pI;
        mesh.globalData().syncPointData
        (
            allPI,
            maxEqOp<scalar>(),
            mapDistribute::transform()
        );
        pI = scalarField(allPI, patch_.meshPoints());
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p
)
:
    patch_(p),
    intact_(nullptr),
    burst_(new burstModels::none())
{}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    intact_(nullptr),
    burst_(burstModel::New(dict))
{}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb
)
:
    patch_(p),
    intact_(nullptr),
    burst_(bppb.burst_->clone())
{}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb,
    const label newSize,
    const label newStart
)
:
    patch_(p),
    intact_(nullptr),
    burst_(bppb.burst_->clone())
{}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb,
    const labelUList& mapAddressing
)
:
    patch_(p),
    intact_(nullptr),
    burst_(bppb.burst_->clone())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPolyPatchBase::~burstPolyPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::burstPolyPatchBase::update
(
    const fvPatch& p,
    scalarField& intact
) const
{
    return returnReduce(burst_->update(p, intact), orOp<bool>());
}


const Foam::Field<Foam::scalar>&
Foam::burstPolyPatchBase::pointIntact() const
{
    if (!pointIntact_.valid())
    {
        makePointIntact();
    }
    return pointIntact_();
}


void Foam::burstPolyPatchBase::write(Ostream& os) const
{
    burst_->writeData(os);
}


// ************************************************************************* //
