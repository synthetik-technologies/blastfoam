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
    usePressure_(false),
    pName_("p"),
    pBurst_(great),
    pRef_(0.0),
    useImpulse_(false),
    impulseName_("impulse"),
    impulseBurst_(great),
    partialBurst_(false),
    curTimeIndex_(-1)
{}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    intact_(nullptr),
    usePressure_(dict.lookupOrDefault("usePressure", true)),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    pBurst_(great),
    pRef_(0.0),
    useImpulse_(dict.lookupOrDefault("useImpulse", false)),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    impulseBurst_(great),
    partialBurst_(dict.lookupOrDefault("partialBurst", false)),
    curTimeIndex_(-1)
{
    if (usePressure_)
    {
        pBurst_ = dict.lookup<scalar>("pBurst");
        if (dict.found("pRef"))
        {
            pRef_ = dict.lookup<scalar>("pRef");
        }
    }
    if (useImpulse_)
    {
        impulseBurst_ = dict.lookup<scalar>("impulseBurst");
    }
}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb
)
:
    patch_(p),
    intact_(nullptr),
    usePressure_(bppb.usePressure_),
    pName_(bppb.pName_),
    pBurst_(bppb.pBurst_),
    pRef_(bppb.pRef_),
    useImpulse_(bppb.useImpulse_),
    impulseName_(bppb.impulseName_),
    impulseBurst_(bppb.impulseBurst_),
    partialBurst_(bppb.partialBurst_),
    curTimeIndex_(-1)
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
    usePressure_(bppb.usePressure_),
    pName_(bppb.pName_),
    pBurst_(bppb.pBurst_),
    pRef_(bppb.pRef_),
    useImpulse_(bppb.useImpulse_),
    impulseName_(bppb.impulseName_),
    impulseBurst_(bppb.impulseBurst_),
    partialBurst_(bppb.partialBurst_),
    curTimeIndex_(-1)
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
    usePressure_(bppb.usePressure_),
    pName_(bppb.pName_),
    pBurst_(bppb.pBurst_),
    pRef_(bppb.pRef_),
    useImpulse_(bppb.useImpulse_),
    impulseName_(bppb.impulseName_),
    impulseBurst_(bppb.impulseBurst_),
    partialBurst_(bppb.partialBurst_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPolyPatchBase::~burstPolyPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::burstPolyPatchBase::update
(
    const scalarField& p,
    const scalarField& impulse,
    scalarField& intact
) const
{
    bool burst = false;
    if (partialBurst_)
    {
        scalarField burstP(intact.size(), -great);
        scalarField burstImp(intact.size(), -great);
        if (usePressure_)
        {
            burstP = p - pBurst_;
        }
        if (useImpulse_)
        {
            burstImp = impulse - impulseBurst_;
        }
        forAll(intact, facei)
        {
            if (intact[facei] > 0.5)
            {
                intact[facei] = burstP[facei] < 0 && burstImp[facei] < 0;
                if (!intact[facei])
                {
                    burst = true;
                }
            }

        }
    }
    else
    {
        // Patch has already burst
        if (gMin(intact) > 0)
        {
            if (usePressure_)
            {
                burst = burst || gMax(p) > pBurst_;
            }
            if (useImpulse_)
            {
                burst = burst || gMax(impulse) > impulseBurst_;
            }
            intact = !burst;
        }
    }
    return returnReduce(burst, orOp<bool>());
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
    writeEntry<Switch>(os, "partialBurst", partialBurst_);

    writeEntry<Switch>(os, "usePressure", usePressure_);
    if (usePressure_)
    {
        writeEntryIfDifferent<word>(os, "pName", "p", pName_);
        writeEntry(os, "pBurst", pBurst_);
        writeEntry(os, "pRef", pRef_);
    }

    writeEntry(os, "useImpulse", useImpulse_);
    if (useImpulse_)
    {
        writeEntryIfDifferent<word>
        (
            os,
            "impulseName",
            "impulse",
            impulseName_
        );
        writeEntry(os, "impulseBurst", impulseBurst_);
    }
}


// ************************************************************************* //
