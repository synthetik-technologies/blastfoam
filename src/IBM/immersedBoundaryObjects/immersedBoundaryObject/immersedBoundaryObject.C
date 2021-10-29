/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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

#include "immersedBoundaryObject.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "immersedBoundaryFvPatchFields.H"
#include "meshSizeObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryObject, 0);
    defineRunTimeSelectionTable(immersedBoundaryObject, dictionary);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::immersedBoundaryObject::interpolationMethod,
        1
    >::names[] =
    {
        "inverseDistance"
    };
}


const
Foam::NamedEnum<Foam::immersedBoundaryObject::interpolationMethod, 1>
    Foam::immersedBoundaryObject::interpolationMethodNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::addIfUnique
(
    List<label>& l,
    const label li,
    const label entry

)
{
    if (li >= l.size())
    {
        l.resize(2*li);
    }
    for (label i = 0; i < li; i++)
    {
        if (l[i] == entry)
        {
            return false;
        }
    }
    l[li] = entry;
    return true;
}


bool Foam::addIfUniqueAndIncrement
(
    List<label>& l,
    label& li,
    const label entry

)
{
    if (li >= l.size())
    {
        l.resize(2*li);
    }
    for (label i = 0; i < li; i++)
    {
        if (l[i] == entry)
        {
            return false;
        }
    }
    l[li++] = entry;
    return true;
}


Foam::scalar Foam::immersedBoundaryObject::delta
(
    const point& p,
    const point& c,
    const scalar& h
) const
{
    if (h < small)
    {
        return 0.0;
    }
    scalar w = 1.0;
    vector d = c - p;
    vector R = shape_->zeroDir(d)/h;
    forAll(p, i)
    {
        scalar r = mag(R[i]);
        if (r <= 0.5)
        {
            w *= 1.0/3.0*(1.0 + sqrt(1.0 - 3.0*sqr(r)));
        }
        else if (r <= 1.5)
        {
            w *= 1.0/6.0*(5.0 - 3.0*r- sqrt(1.0 - 3.0*sqr(1.0 - r)));
        }
        else
        {
            return 0.0;
        }
    }
    return w;
}


void Foam::immersedBoundaryObject::calcMapping() const
{
    const pointField& points = pMesh_.points();

    interpToCells_.resize(nFaces());
    interpToWeights_.resize(nFaces());
    deltaCoeffs_.resize(nFaces());

    boundBox meshBb(points, false); // local
    meshBb.inflate(1e-4);
    bool notFound =
        !shape_->bounds().overlaps(meshBb)
     && !meshBb.contains(shape_->bounds())
     && !shape_->bounds().contains(meshBb);
    if (notFound)
    {
        internalCellsPtr_ = new labelList();
        boundaryCellsPtr_ = new labelList();
        shellCellsPtr_ = new labelList();

        interpFromPoints_.clear();
        interpFromWeights_.clear();
        internalC_.clear();

        patchInternalCells_ = -1;
        patchExternalCells_ = -1;

        interpToCells_ = labelList();
        interpToWeights_ = scalarList();
        deltaCoeffs_ = 0.0;
    }

    if (returnReduce(notFound, andOp<bool>()))
    {
        return;
    }

    const volScalarField& edgeLength = meshSizeObject::New(pMesh_).dx();

    labelList insidePoints(shape_->calcInside(points));
    boolList inside(points.size(), false);
    forAll(insidePoints, i)
    {
        inside[insidePoints[i]] = true;
    }

    internalCellsPtr_ = new labelList(pMesh_.nCells());
    boundaryCellsPtr_ = new labelList(pMesh_.nCells());


    labelList& internalCells = *internalCellsPtr_;
    labelList& boundaryCells = *boundaryCellsPtr_;

    PackedBoolList pbShellCells(pMesh_.nCells(), false);

    label ii = 0;
    label bi = 0;
    forAll(pMesh_.cells(), celli)
    {
        const labelList& cp = pMesh_.cellPoints()[celli];
        bool allIn = true;
        bool allOut = true;
        forAll(cp, i)
        {
            allIn = allIn && inside[cp[i]];
            allOut = allOut && !inside[cp[i]];
        }
        if (allIn)
        {
            internalCells[ii++] = celli;
        }
        else if (!allIn && !allOut)
        {
            boundaryCells[bi++] = celli;
            pbShellCells.set(celli, true);
        }
    }
    internalCells.resize(ii);
    boundaryCells.resize(bi);

    PackedBoolList pbShellCellsOrig(pbShellCells);
    forAll(pbShellCellsOrig, celli)
    {
        if (pbShellCellsOrig.get(celli))
        {
            const List<label>& neighbours =
                cellNeighbours_.cellCell(celli).localStencil();
            forAll(neighbours, cellj)
            {
                pbShellCells.set(neighbours[cellj], true);
            }
        }
    }

    label si = 0;
    shellCellsPtr_ = new labelList(pbShellCells.count());

    labelList& shellCells(*shellCellsPtr_);
    forAll(pbShellCells, i)
    {
        if (pbShellCells.get(i))
        {
            shellCells[si++] = i;
        }
    }

    labelHashSet insideCells;
    labelHashSet setCells;

    patchInternalCells_ = labelList(nFaces(), -1);
    patchExternalCells_ = labelList(nFaces(), -1);
    Map<label> shellCellOwners;

    const globalIndex& gIndex(cellNeighbours_.globalNumbering());

    //- Find global owner of the face centre
    labelList ownerCell(nFaces(), -1);
    scalarList eL(nFaces(), great);
    forAll(patch_, facei)
    {
        const point& fc = patch_.faceCentres()[facei];
        label celli = -1;
//         if (interpToCells_[facei].size())
//         {
//             if (interpToCells_[facei][0] >= 0)
//             {
//                 if
//                 (
//                     pMesh_.primitiveMesh::pointInCell(fc, interpToCells_[facei][0])
//                 )
//                 {
//                     celli = interpToCells_[facei][0];
//                 }
//             }
//         }
        if (meshBb.contains(fc))// && celli < 0)
        {
            celli = pMesh_.findCell(fc, polyMesh::FACE_PLANES);
        }
        if (celli != -1)
        {
            eL[facei] = edgeLength[celli];
            ownerCell[facei] = gIndex.toGlobal(celli);
        }
    }
    reduce(ownerCell, maxOp<labelList>());
    reduce(eL, minOp<scalarList>());

    const vectorField& Sf(this->Sf());
    scalarField magSf(this->magSf());

    scalarList maxIntNDotN(nFaces(), 0.0);
    scalarList maxExtNDotN(nFaces(), 0.0);
    labelList intCell(ownerCell);
    labelList extCell(ownerCell);
    scalarList sumWs(nFaces(), 0.0);
    forAll(patch_, facei)
    {
        const point& fc = patch_.faceCentres()[facei];
        const label gCelli = ownerCell[facei];
        if (!cellNeighbours_.cellCellMap().found(ownerCell[facei]) || gCelli == -1)
        {
            interpToCells_[facei].clear();
            interpToWeights_[facei].clear();
            patchInternalCells_[facei] = -1;
            patchExternalCells_[facei] = -1;
            deltaCoeffs_[facei] = 0.0;
            continue;
        }
        const cellStencil& stencil
        (
            cellNeighbours_.cellCellMap()[gCelli]
        );
        const vector& cci = stencil.centre();
        const List<label>& neighbors = stencil.localStencil();

        label celli = gIndex.isLocal(gCelli) ? gIndex.toLocal(gCelli) : -1;
        List<label> interpToCells(neighbors.size(), -1);
        List<scalar> interpToWeights(neighbors.size(), 0.0);
        label fi = 0;
        if (celli != -1)
        {
            interpToCells[fi] = celli;
            interpToWeights[fi] = delta(fc, cci, eL[facei]);
            shellCellOwners.insert(celli, facei);
            fi++;
        }

        const bool in(shape_->inside(cci));
        setCells.insert(gCelli);
        vector n(Sf[facei]/magSf[facei]);
        vector diff(fc - cci);
        diff /= mag(diff);
        scalar nDotN(mag(diff & n));
        if (in)
        {
            insideCells.set(gCelli);
            intCell[facei] = celli;
            extCell[facei] = -1;
            maxIntNDotN[facei] = nDotN;
        }
        else
        {
            intCell[facei] = -1;
            extCell[facei] = celli;
            maxExtNDotN[facei] = nDotN;
        }

        //- If the cell is not inside the object set as the internal cells
        forAll(neighbors, j)
        {
            const label cellj = neighbors[j];
            const label gCellj = gIndex.toGlobal(cellj);
            const vector& ccj = pMesh_.cellCentres()[cellj];

            if (addIfUnique(interpToCells, fi, cellj))
            {
                scalar w = delta(fc, ccj, eL[facei]);
                interpToWeights[fi] = w;
                fi++;
            }
            if (!shellCellOwners.found(cellj))
            {
                shellCellOwners.insert(cellj, facei);
            }

            bool inj;
            if (!setCells.found(gCelli))
            {
                inj = shape_->inside(ccj);
                setCells.insert(gCellj);
            }
            else
            {
                inj = insideCells.found(gCellj);
            }

            diff = fc - ccj;
            diff /= mag(diff);
            nDotN = mag(diff & n);
            if (inj && !in)
            {
                insideCells.insert(gCellj);
                if (nDotN > maxIntNDotN[facei])
                {
                    intCell[facei] = cellj;
                    maxIntNDotN[facei] = nDotN;
                }
            }
            else if (!inj && in)
            {
                if (nDotN > maxExtNDotN[facei])
                {
                    extCell[facei] = cellj;
                    maxExtNDotN[facei] = nDotN;
                }
            }

        }
        interpToCells.resize(fi);
        interpToWeights.resize(fi);
        sumWs[facei] = sum(interpToWeights);
        interpToCells_[facei].transfer(interpToCells);
        interpToWeights_[facei].transfer(interpToWeights);
    }


    scalarList maxIntD(returnReduce(maxIntNDotN, maxOp<scalarList>()));
    scalarList maxExtD(returnReduce(maxExtNDotN, maxOp<scalarList>()));
    forAll(patch_, facei)
    {
        if (maxIntD[facei] == maxIntNDotN[facei])
        {
            patchInternalCells_[facei] = intCell[facei];
        }
        else
        {
            patchInternalCells_[facei] = -1;
        }

        if (maxExtD[facei] == maxExtNDotN[facei])
        {
            patchExternalCells_[facei] = extCell[facei];
        }
        else
        {
            patchExternalCells_[facei] = -1;
        }

        if (patchExternalCells_[facei] >= 0)
        {
            deltaCoeffs_[facei] =
                1.0
                /mag
                (
                    shape().zeroDir
                    (
                        pMesh_.cellCentres()[patchExternalCells_[facei]]
                      - patch_.faceCentres()[facei]
                    )
                );
        }
        else
        {
            deltaCoeffs_[facei] = 0.0;
        }
    }

    reduce(sumWs, sumOp<scalarList>());
    forAll(patch_, facei)
    {
        forAll(interpToWeights_[facei], j)
        {
            interpToWeights_[facei][j] /= sumWs[facei];
        }
    }

    interpFromPoints_.resize(shellCellsPtr_->size());
    interpFromWeights_.resize(shellCellsPtr_->size());

    forAll(shellCells, i)
    {
        labelList interpFromPoints(nFaces(), -1);
        scalarList interpFromWeights(nFaces(), 0.0);
        const label celli = shellCells[i];
        const vector& x = pMesh_.cellCentres()[celli];
        label nPts = 0;
        scalar sumW = 0;

        if (!shellCellOwners.found(celli))
        {
            continue;
        }
        const labelList& nearest = nearestNeighbours_[shellCellOwners[celli]];
        forAll(nearest, fi)
        {
            label facei = nearest[fi];
            if (!interpToCells_[facei].size())
            {
                continue;
            }

            const vector& fc = patch_.faceCentres()[facei];
            if (mag(fc - x) < 1.5*edgeLength[celli])
            {
                interpFromPoints[nPts] = facei;
                scalar d = delta(x, fc, eL[facei]);
                scalar invD(1.0/max(mag(x - fc), small));
                interpFromWeights[nPts] = d*invD;
                sumW += invD;
                nPts++;
            }
        }
        interpFromPoints.resize(nPts);
        interpFromWeights.resize(nPts);
        forAll(interpFromWeights, j)
        {
            interpFromWeights[j] /= max(sumW, small);
        }
        interpFromPoints_[i].transfer(interpFromPoints);
        interpFromWeights_[i].transfer(interpFromWeights);
    }

    internalC_.resize(internalCells.size());
    forAll(internalC_, i)
    {
        internalC_[i] = pMesh_.cellCentres()[internalCells[i]];
    }
}

void Foam::immersedBoundaryObject::updateWeights() const
{
    const volScalarField& edgeLength = meshSizeObject::New(pMesh_).dx();
    forAll(patch_, facei)
    {
        const point& fc = patch_.faceCentres()[facei];

        if (!interpToCells_[facei].size())
        {
            continue;
        }
        label celli = interpToCells_[facei][0];

        scalar sumW = 0.0;

        const labelList& neighbours = interpToCells_[facei];
        forAll(neighbours, j)
        {
            const label cellj = neighbours[j];
            scalar w =
                delta
                (
                    fc,
                    pMesh_.cellCentres()[cellj],
                    edgeLength[celli]
                );
            interpToWeights_[facei][j] = w;
            sumW += w;
        }
        forAll(interpToWeights_[facei], j)
        {
            interpToWeights_[facei][j] /= max(sumW, small);
        }
        deltaCoeffs_[facei] =
            1.0
            /mag
            (
                shape().zeroDir
                (
                    pMesh_.cellCentres()[patchExternalCells_[facei]] - fc
                )
            );
    }

    forAll(interpFromPoints_, i)
    {
        const label celli = (*shellCellsPtr_)[i];
        const vector& x = pMesh_.cellCentres()[celli];

        const labelList& nearest = interpFromPoints_[i];
        scalar sumW = 0.0;
        forAll(nearest, fi)
        {
            const label own = interpToCells_[fi][0];
            const vector& fc = patch_.faceCentres()[fi];
            if (mag(fc - x) < 1.5*edgeLength[celli])
            {
                scalar d = delta(x, fc, edgeLength[own]);
                scalar invD(1.0/max(mag(x - fc), small));
                interpFromWeights_[i][fi] = d*invD;
                sumW += invD;
            }
        }
        forAll(interpFromWeights_[i], j)
        {
            interpFromWeights_[i][j] /= max(sumW, small);
        }
    }
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObject::immersedBoundaryObject
(
    const polyMesh& pMesh,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    pMesh_(pMesh),
    cellNeighbours_(extendedNLevelCPCCellToCellStencil::New(pMesh)),
    dict_(dict),
    geometricD_(1, 1, 1),
    shape_(immersedShape::New(pMesh, *this, dict)),
    mode_(INVERSE_DISTANCE),
    patch_(shape_->patch()),
    boundaryCellsPtr_(nullptr),
    internalCellsPtr_(nullptr),
    shellCellsPtr_(nullptr),
    patchInternalCells_(nFaces(), -1),
    patchExternalCells_(nFaces(), -1),
    internalC_(),
    interpToCells_(nFaces()),
    interpToWeights_(nFaces()),
    interpFromPoints_(nFaces()),
    interpFromWeights_(nFaces()),
    deltaCoeffs_(nFaces()),
    nearestNeighbours_(nFaces()),
    mass_(great),
    force_(nFaces()),
    forceEff_(Zero),
    momentEff_(Zero),
    forceExt_(Zero),
    momentExt_(Zero),
    g_
    (
        pMesh_.template foundObject<uniformDimensionedVectorField>("g")
      ? pMesh_.template lookupObject<uniformDimensionedVectorField>("g").value()
      : Zero
    ),
    report_(dict.lookupOrDefault<Switch>("report", true)),
    scalarBoundaries_(0),
    vectorBoundaries_(0)
{
    // Find largest grid spacing
    scalar maxEdgeLength =
        max(meshSizeObject::New(pMesh).dx()).value();

    // Scale by 1.5 for weighting
    maxEdgeLength *= 1.5;

    // Add face centers within the maxEdgeLength to the map
    const pointField& faceCentres(patch_.faceCentres());
    forAll(nearestNeighbours_, fi)
    {
        List<label> nn(10, -1);
        label ni = 0;
        scalarList diff(mag(faceCentres[fi] - faceCentres));
        forAll(diff, fj)
        {
            if (ni == nn.size())
            {
                nn.resize(nn.size()*2);
            }
            if (diff[fj] < maxEdgeLength)
            {
                nn[ni++] = fj;
            }
        }
        nn.resize(ni);
        nearestNeighbours_[fi].transfer(nn);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObject::~immersedBoundaryObject()
{
    clearOut();
}


void Foam::immersedBoundaryObject::clearOut()
{
    deleteDemandDrivenData(boundaryCellsPtr_);
    deleteDemandDrivenData(internalCellsPtr_);
    deleteDemandDrivenData(shellCellsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedBoundaryObject::setValues()
{
    forAll(scalarBoundaries_, i)
    {
        scalarBoundaries_[i].setValues();
    }
    forAll(vectorBoundaries_, i)
    {
        vectorBoundaries_[i].setValues();
    }
}


template<>
void Foam::immersedBoundaryObject::addField
(
    volScalarField& f,
    const dictionary& defaultDict
)
{
    dictionary dict(f.name());
    bool found = false;
    if (this->dict_.found("boundaries"))
    {
        if (this->dict_.subDict("boundaries").found(f.name()))
        {
            dict = this->dict_.subDict("boundaries").subDict(f.name());
            found = true;
        }
    }
    if (!found)
    {
        dict = defaultDict.subDict(f.name());
    }
    scalarBoundaries_.resize(scalarBoundaries_.size() + 1);
    scalarBoundaries_.set
    (
        scalarBoundaries_.size() - 1,
        f.name(),
        immersedBoundaryScalarPatchField::New(f, dict, *this).ptr()
    );
}


template<>
void Foam::immersedBoundaryObject::addForcing
(
    const word& name,
    volScalarField& F,
    const volScalarField& alphaRho,
    const volScalarField& alphaRhoFOld,
    const volScalarField& RHS,
    const dimensionedScalar& dt
) const
{
    scalarBoundaries_[name].addForcing
    (
        F,
        alphaRho,
        alphaRhoFOld,
        RHS,
        dt.value()
    );
}

template<>
void Foam::immersedBoundaryObject::addField
(
    volVectorField& f,
    const dictionary& defaultDict
)
{
    dictionary dict(f.name());
    bool found = false;
    if (this->dict_.found("boundaries"))
    {
        if (this->dict_.subDict("boundaries").found(f.name()))
        {
            dict = this->dict_.subDict("boundaries").subDict(f.name());
            found = true;
        }
    }
    if (!found)
    {
        dict = defaultDict.subDict(f.name());
    }
    vectorBoundaries_.resize(vectorBoundaries_.size() + 1);
    vectorBoundaries_.set
    (
        vectorBoundaries_.size() - 1,
        f.name(),
        immersedBoundaryVectorPatchField::New(f, dict, *this).ptr()
    );
}


template<>
void Foam::immersedBoundaryObject::addForcing
(
    const word& name,
    volVectorField& F,
    const volScalarField& alphaRho,
    const volVectorField& alphaRhoFOld,
    const volVectorField& RHS,
    const dimensionedScalar& dt
) const
{
    vectorBoundaries_[name].addForcing
    (
        F,
        alphaRho,
        alphaRhoFOld,
        RHS,
        dt.value()
    );
}


const Foam::fvMesh& Foam::immersedBoundaryObject::immersedMesh() const
{
    NotImplemented;
    return dynamic_cast<const fvMesh&>(pMesh_);
}


const Foam::immersedMeshMapper* Foam::immersedBoundaryObject::mapper() const
{
    NotImplemented;
    return nullptr;
}


void Foam::immersedBoundaryObject::status() const
{
    Info<< name() << ":" <<endl;
}


Foam::labelList Foam::immersedBoundaryObject::calcInside
(
    const pointField& points
) const
{
    return shape_->calcInside(points);
}



Foam::tmp<Foam::vectorField>
Foam::immersedBoundaryObject::velocity() const
{
    return
        (faceCentres() - faceCentresOld())
       /pMesh_.time().deltaTValue();
}


Foam::vector
Foam::immersedBoundaryObject::velocity(const label facei) const
{
    return
        (faceCentres()[facei] - faceCentresOld()[facei])
       /pMesh_.time().deltaTValue();
}


bool Foam::immersedBoundaryObject::read(const dictionary& dict)
{
    report_ = dict.lookupOrDefault<Switch>("report", false);
    return true;
}


void Foam::immersedBoundaryObject::write(Ostream& os) const
{
    writeEntry(os, "report", report_);
}

void Foam::immersedBoundaryObject::write(dictionary& dict) const
{
    dict.add("report", report_);
}

// ************************************************************************* //
