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
#include "uniformDimensionedFields.H"
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
    if (!shape_.valid())
    {
        initialize();
    }
    const pointField& points = pMesh_.points();

    interpToCells_.resize(this->size());
    interpToWeights_.resize(this->size());
    deltaCoeffs_.resize(this->size());

    boundBox meshBb(points, false); // local
    meshBb.inflate(1e-4);
    bool notFound =
        !shape_->bounds().overlaps(meshBb)
     && !meshBb.contains(shape_->bounds())
     && !shape_->bounds().contains(meshBb);
    if (notFound)
    {
        allInternalCellsPtr_ = new labelList();
        internalCellsPtr_ = new labelList();
        boundaryCellsPtr_ = new labelList();
        shellCellsPtr_ = new labelList();

        interpFromPoints_.clear();
        interpFromWeights_.clear();

        patchInternalCells_ = -1;
        patchExternalCells_ = -1;
        patchCells_.clear();
        patchMap_.clear();

        interpToCells_ = labelList();
        interpToWeights_ = scalarList();
        deltaCoeffs_ = great;
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
    boundaryCellsPtr_ = new labelList(pMesh_.nCells());


    labelHashSet internalCells;
    labelList& boundaryCells = *boundaryCellsPtr_;

    PackedBoolList pbShellCells(pMesh_.nCells(), false);

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
            internalCells.insert(celli);
        }
        else if (!allIn && !allOut)
        {
            boundaryCells[bi++] = celli;
            pbShellCells.set(celli, true);
        }
    }
    allInternalCellsPtr_ = new labelList(internalCells.toc());
    boundaryCells.resize(bi);

    const extendedNLevelCPCCellToCellStencil& cellNeighbours =
        extendedNLevelCPCCellToCellStencil::New(pMesh_);

    PackedBoolList pbShellCellsOrig(pbShellCells);
    forAll(pbShellCellsOrig, celli)
    {
        if (pbShellCellsOrig.get(celli))
        {
            const List<label>& neighbours =
                cellNeighbours.cellCell(celli).localStencil(2);
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

    patchInternalCells_ = labelList(this->size(), -1);
    patchExternalCells_ = labelList(this->size(), -1);
    Map<label> shellCellOwners;

    const globalIndex& gIndex(cellNeighbours.globalNumbering());

    //- Find global owner of the face centre
    labelList ownerCell(this->size(), -1);
    scalarList eL(this->size(), great);
    const vectorField& faceCentres(this->faceCentres());
    forAll(*this, facei)
    {
        const point& fc = faceCentres[facei];
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

    const vectorField& Sf(this->faceAreas());
    scalarField magSf(mag(Sf));

    scalarList minInt(this->size(), great);
    scalarList minExt(this->size(), great);
    labelList intCell(ownerCell.size(), -1);
    labelList extCell(ownerCell.size(), -1);
    scalarList sumWs(this->size(), 0.0);
    forAll(*this, facei)
    {
        const point& fc = faceCentres[facei];
        const label gCelli = ownerCell[facei];
        if
        (
            !cellNeighbours.cellCellMap().found(ownerCell[facei])
         || gCelli == -1
        )
        {
            interpToCells_[facei].clear();
            interpToWeights_[facei].clear();
            patchInternalCells_[facei] = -1;
            patchExternalCells_[facei] = -1;
            deltaCoeffs_[facei] = great;
            continue;
        }
        const cellStencil& stencil
        (
            cellNeighbours.cellCellMap()[gCelli]
        );
        const vector& cci = stencil.centre();
        const List<label> neighbors(stencil.localStencil(0, 2));

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
        vector diff(fc - cci);
        if (in)
        {
            insideCells.set(gCelli);
            intCell[facei] = celli;
            minInt[facei] = mag(diff);
        }
        else
        {
            extCell[facei] = celli;
            minExt[facei] = mag(diff);
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
            if (!setCells.found(gCellj))
            {
                inj = shape_->inside(ccj);
                setCells.insert(gCellj);
                if (inj)
                {
                    insideCells.insert(gCellj);
                }
            }
            else
            {
                inj = insideCells.found(gCellj);
            }

            diff = fc - ccj;
            scalar magDiff(mag(diff));
            if (inj && magDiff < minInt[facei])
            {
                intCell[facei] = cellj;
                minInt[facei] = magDiff;
            }
            else if (!inj && magDiff < minExt[facei])
            {
                extCell[facei] = cellj;
                minExt[facei] = magDiff;
            }

        }
        interpToCells.resize(fi);
        interpToWeights.resize(fi);
        sumWs[facei] = sum(interpToWeights);
        interpToCells_[facei].transfer(interpToCells);
        interpToWeights_[facei].transfer(interpToWeights);
    }

    scalarList minIntG(returnReduce(minInt, minOp<scalarList>()));
    scalarList minExtG(returnReduce(minExt, minOp<scalarList>()));

    Map<DynamicList<label>> patchMap;
    forAll(*this, facei)
    {
        if (minIntG[facei] == minInt[facei])
        {
            patchInternalCells_[facei] = intCell[facei];
            internalCells.erase(intCell[facei]);
        }
        else
        {
            patchInternalCells_[facei] = -1;
        }

        if (minExtG[facei] == minExt[facei])
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
                      - this->faceCentres()[facei]
                    )
                );
        }
        else
        {
            deltaCoeffs_[facei] = great;
        }

        label celli = -1;
        if (patchExternalCells_[facei] >= 0)
        {
            celli = patchExternalCells_[facei];
        }
        else if (patchInternalCells_[facei] >= 0)
        {
            celli = patchInternalCells_[facei];
        }
        if (celli >= 0)
        {
            if (!patchMap.found(celli))
            {
                patchMap.insert(celli, DynamicList<label>(1, facei));
            }
            else
            {
                patchMap[celli].append(facei);
            }
        }
    }
    patchCells_.resize(patchMap.size());
    patchMap_.resize(patchMap.size());
    label pci = 0;
    forAllIter(Map<DynamicList<label>>, patchMap, iter)
    {
        patchCells_[pci] = iter.key();
        patchMap_[pci++].transfer(iter());
    }

    reduce(sumWs, sumOp<scalarList>());
    forAll(*this, facei)
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
        labelList interpFromPoints(this->size(), -1);
        scalarField interpFromWeights(this->size(), 0.0);
        const label celli = shellCells[i];
        const vector& x = pMesh_.cellCentres()[celli];
        label nPts = 0;
        scalar sumW = 0;

        if (!shellCellOwners.found(celli))
        {
            interpFromPoints_[i].clear();
            interpFromWeights_[i].clear();
            continue;
        }
        const labelList& nearest = nearestNeighbours_[shellCellOwners[celli]];
        forAll(nearest, fi)
        {
            label facei = nearest[fi];

            const vector& fc = this->faceCentres()[facei];
            if (mag(fc - x) < 1.5*eL[facei])
            {
                interpFromPoints[nPts] = facei;
                scalar d = delta(x, fc, eL[facei]);
                interpFromWeights[nPts] = d;
                sumW += d;
                nPts++;
            }
        }
        interpFromPoints.resize(nPts);
        interpFromWeights.resize(nPts);
        interpFromWeights /= max(sumW, small);

        interpFromPoints_[i].transfer(interpFromPoints);
        interpFromWeights_[i].transfer(interpFromWeights);
    }

    internalCellsPtr_ = new labelList(internalCells.toc());
}


// void Foam::immersedBoundaryObject::updateWeights() const
// {
//     const volScalarField& edgeLength = meshSizeObject::New(pMesh_).dx();
//     forAll(*this, facei)
//     {
//         const point& fc = this->faceCentres()[facei];
//
//         if (!interpToCells_[facei].size())
//         {
//             continue;
//         }
//         label celli = interpToCells_[facei][0];
//
//         scalar sumW = 0.0;
//
//         const labelList& neighbours = interpToCells_[facei];
//         forAll(neighbours, j)
//         {
//             const label cellj = neighbours[j];
//             scalar w =
//                 delta
//                 (
//                     fc,
//                     pMesh_.cellCentres()[cellj],
//                     edgeLength[celli]
//                 );
//             interpToWeights_[facei][j] = w;
//             sumW += w;
//         }
//         forAll(interpToWeights_[facei], j)
//         {
//             interpToWeights_[facei][j] /= max(sumW, small);
//         }
//         deltaCoeffs_[facei] =
//             1.0
//             /mag
//             (
//                 shape().zeroDir
//                 (
//                     pMesh_.cellCentres()[patchExternalCells_[facei]] - fc
//                 )
//             );
//     }
//
//     forAll(interpFromPoints_, i)
//     {
//         const label celli = (*shellCellsPtr_)[i];
//         const vector& x = pMesh_.cellCentres()[celli];
//
//         const labelList& nearest = interpFromPoints_[i];
//         scalar sumW = 0.0;
//         forAll(nearest, fi)
//         {
//             const label own = interpToCells_[fi][0];
//             const vector& fc = this->faceCentres()[fi];
//             if (mag(fc - x) < 1.5*edgeLength[celli])
//             {
//                 scalar d = delta(x, fc, edgeLength[own]);
//                 scalar invD(1.0/max(mag(x - fc), small));
//                 interpFromWeights_[i][fi] = d*invD;
//                 sumW += d*invD;
//             }
//         }
//         forAll(interpFromWeights_[i], j)
//         {
//             interpFromWeights_[i][j] /= max(sumW, small);
//         }
//     }
// }


void Foam::immersedBoundaryObject::initialize() const
{
    //- Set the stencil to be 2 levels deep
    extendedNLevelCPCCellToCellStencil::New(pMesh_).setNLevel(2);

    shape_ = immersedShape::New(pMesh_, *this, dict_);

    const standAlonePatch& csap(*this);
    standAlonePatch& sap(const_cast<standAlonePatch&>(csap));

    sap = move(shape_->createPatch()());

    sap.movePoints(this->transform(this->points()));
    shape_->movePoints(this->points());

    patchInternalCells_.resize(this->size(), -1);
    patchExternalCells_.resize(this->size(), -1);
    interpToCells_.resize(this->size());
    interpToWeights_.resize(this->size());
    interpFromPoints_.resize(this->size());
    interpFromWeights_.resize(this->size());
    deltaCoeffs_.resize(this->size());
    nearestNeighbours_.resize(this->size());
    force_.setSize(this->size());
    force_ = Zero;

    // Find largest grid spacing
    scalar maxEdgeLength =
        max(meshSizeObject::New(pMesh_).dx()).value();

    // Scale by 1.5 for weighting
    maxEdgeLength *= 1.5;

    // Add face centers within the maxEdgeLength to the map
    const pointField& faceCentres(sap.faceCentres());
    forAll(nearestNeighbours_, fi)
    {
        DynamicList<label> nn;
        scalarList diff(mag(faceCentres[fi] - faceCentres));
        forAll(diff, fj)
        {
            if (diff[fj] < maxEdgeLength)
            {
                nn.append(fj);
            }
        }
        nearestNeighbours_[fi].transfer(nn);
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObject::immersedBoundaryObject
(
    const polyPatch& patch,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    standAlonePatch(faceList(), pointField()),
    patch_(patch),
    pMesh_(patch.boundaryMesh().mesh()),
    dict_(dict),
    mode_(INVERSE_DISTANCE),
    shape_(nullptr),
    boundaryCellsPtr_(nullptr),
    allInternalCellsPtr_(nullptr),
    internalCellsPtr_(nullptr),
    shellCellsPtr_(nullptr),
    patchInternalCells_(),
    patchExternalCells_(),
    interpToCells_(),
    interpToWeights_(),
    interpFromPoints_(),
    interpFromWeights_(),
    deltaCoeffs_(),
    magSfPtr_(nullptr),
    nearestNeighbours_(),
    force_(),
    forceEff_(Zero),
    momentEff_(Zero),
    forceExt_(Zero),
    momentExt_(Zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObject::~immersedBoundaryObject()
{
    clearOut();
}


void Foam::immersedBoundaryObject::clearOut()
{
    deleteDemandDrivenData(boundaryCellsPtr_);
    deleteDemandDrivenData(allInternalCellsPtr_);
    deleteDemandDrivenData(internalCellsPtr_);
    deleteDemandDrivenData(shellCellsPtr_);
    deleteDemandDrivenData(magSfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::immersedBoundaryObject::status(const bool print) const
{
    if (!print)
    {
        return;
    }
    Info<< name() << ":" <<endl;
}


Foam::labelList Foam::immersedBoundaryObject::calcInside
(
    const pointField& points
) const
{
    return shape_->calcInside(points);
}


void Foam::immersedBoundaryObject::movePoints()
{
    standAlonePatch::movePoints(this->transform(this->points()));
    shape_->movePoints(this->points());
    clearOut();
}


bool Foam::immersedBoundaryObject::read(const dictionary& dict)
{
    return true;
}


void Foam::immersedBoundaryObject::write(Ostream& os) const
{}


void Foam::immersedBoundaryObject::write(dictionary& dict) const
{}


// ************************************************************************* //
