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
#include "globalPolyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstPolyPatchDistributor, 0);
    defineTypeNameAndDebug(burstPolyPatchBase, 0);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::burstPolyPatchDistributor::burstPolyPatchDistributor
(
    const polyMesh& mesh
)
:
    BalanceMeshObject(typeName, mesh),
    mesh_(mesh)
{}

Foam::burstPolyPatchDistributor& Foam::burstPolyPatchDistributor::New
(
    const polyMesh& mesh
)
{
    if (!mesh.foundObject<burstPolyPatchDistributor>(typeName))
    {
        burstPolyPatchDistributor* burstList =
            new burstPolyPatchDistributor(mesh);
        burstList->store(burstList);
    }
    return mesh.lookupObjectRef<burstPolyPatchDistributor>(typeName);
}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p
)
:
    patch_(p),
    intact_(p.size(), 1.0),
    usePressure_(false),
    pName_("p"),
    pBurst_(great),
    pRef_(0.0),
    useImpulse_(false),
    impulseName_("impulse"),
    impulseBurst_(great),
    partialBurst_(false),
    needRead_(false),
    curTimeIndex_(-1)
{
    burstPolyPatchDistributor::New(patch_.boundaryMesh().mesh());
}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    intact_(p.size(), 1.0),//"intact", dict, p.size()),
    usePressure_(dict.lookupOrDefault("usePressure", true)),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    pBurst_(great),
    pRef_(0.0),
    useImpulse_(dict.lookupOrDefault("useImpulse", false)),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    impulseBurst_(great),
    partialBurst_(dict.lookupOrDefault("partialBurst", false)),
    needRead_(false),
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
    burstPolyPatchDistributor::New(patch_.boundaryMesh().mesh());
}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb
)
:
    patch_(p),
    intact_(bppb.intact_),
    usePressure_(bppb.usePressure_),
    pName_(bppb.pName_),
    pBurst_(bppb.pBurst_),
    pRef_(bppb.pRef_),
    useImpulse_(bppb.useImpulse_),
    impulseName_(bppb.impulseName_),
    impulseBurst_(bppb.impulseBurst_),
    partialBurst_(bppb.partialBurst_),
    needRead_(false),
    curTimeIndex_(-1)
{
    burstPolyPatchDistributor::New(patch_.boundaryMesh().mesh());
}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb,
    const label newSize,
    const label newStart
)
:
    patch_(p),
    intact_(bppb.intact_),
    usePressure_(bppb.usePressure_),
    pName_(bppb.pName_),
    pBurst_(bppb.pBurst_),
    pRef_(bppb.pRef_),
    useImpulse_(bppb.useImpulse_),
    impulseName_(bppb.impulseName_),
    impulseBurst_(bppb.impulseBurst_),
    partialBurst_(bppb.partialBurst_),
    needRead_(false),
    curTimeIndex_(-1)
{
    burstPolyPatchDistributor::New(patch_.boundaryMesh().mesh());
}


Foam::burstPolyPatchBase::burstPolyPatchBase
(
    const polyPatch& p,
    const burstPolyPatchBase& bppb,
    const labelUList& mapAddressing
)
:
    patch_(p),
    intact_(bppb.intact_, mapAddressing),
    usePressure_(bppb.usePressure_),
    pName_(bppb.pName_),
    pBurst_(bppb.pBurst_),
    pRef_(bppb.pRef_),
    useImpulse_(bppb.useImpulse_),
    impulseName_(bppb.impulseName_),
    impulseBurst_(bppb.impulseBurst_),
    partialBurst_(bppb.partialBurst_),
    needRead_(false),
    curTimeIndex_(-1)
{
    burstPolyPatchDistributor::New(patch_.boundaryMesh().mesh());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPolyPatchBase::~burstPolyPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstPolyPatchDistributor::preDistribute()
{
    intact_.set(new scalarField(mesh_.nFaces(), 1.0));
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<burstPolyPatchBase>(mesh_.boundaryMesh()[patchi]))
        {
            const burstPolyPatchBase& p =
                dynamicCast<const burstPolyPatchBase>
                (
                    mesh_.boundaryMesh()[patchi]
                );
            const label start = p.patch().start();
            forAll(p.patch(), fi)
            {
                const label facei = start + fi;
                intact_()[facei] = p.intact()[fi];
            }
        }
    }
}


void Foam::burstPolyPatchDistributor::distribute
(
    const mapDistributePolyMesh& map
)
{
    Pout<<"distibuting: "<<endl;
    map.distributeFaceData(intact_());
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<burstPolyPatchBase>(mesh_.boundaryMesh()[patchi]))
        {
            burstPolyPatchBase& p =
                const_cast<burstPolyPatchBase&>
                (
                    dynamicCast<const burstPolyPatchBase>
                    (
                        mesh_.boundaryMesh()[patchi]
                    )
                );
            const label start = p.patch().start();
            p.intact_.resize(p.patch().size());
            forAll(p.patch(), fi)
            {
                const label facei = start + fi;
                p.intact_[fi] = intact_()[facei];
            }
        }
    }
    intact_.clear();
}


bool Foam::burstPolyPatchBase::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    bool burst = false;
    if (partialBurst_)
    {
        scalarField burstP(intact_.size(), -great);
        scalarField burstImp(intact_.size(), -great);
        if (usePressure_)
        {
            burstP = p - pBurst_;
        }
        if (useImpulse_)
        {
            burstImp = impulse - impulseBurst_;
        }
        forAll(intact_, facei)
        {
            if (intact_[facei])
            {
                intact_[facei] = burstP[facei] < 0 && burstImp[facei] < 0;
                if (!intact_[facei])
                {
                    burst = true;
                }
            }

        }
    }
    else
    {
        // Patch has already burst
        if (gMin(intact_) > 0)
        {
            if (usePressure_)
            {
                burst = gMax(p) > pBurst_;
            }
            if (useImpulse_)
            {
                burst = gMax(impulse) > impulseBurst_;
            }
            intact_ = !burst;
        }
    }
    return returnReduce(burst, orOp<bool>());
}


const Foam::Field<Foam::scalar>&
Foam::burstPolyPatchBase::pointIntact() const
{
    if (!pointIntact_.valid())
    {
        pointIntact_.set
        (
            globalPolyBoundaryMesh::New
            (
                dynamicCast<const polyMesh>
                (
                    patch_.boundaryMesh().mesh().thisDb()
                )
            )[patch_].faceToPoint(intact_).ptr()
        );
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
