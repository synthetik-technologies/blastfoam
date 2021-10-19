/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Generalized grid interface (GGI) patch, providing coupling
    between arbitrary patches which belong to the same fvMesh

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "ggiFvPatch.H"
#include "fvsPatchFields.H"
#include "fvBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, ggiFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ggiFvPatch::~ggiFvPatch()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Make patch weighting factors
void Foam::ggiFvPatch::makeWeights(scalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (ggiPolyPatch_.owner())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.
        const vectorField n(nf());

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  HJ, 24/Aug/2011
        const scalarField nfc
        (
            mag(n & (ggiPolyPatch_.reconFaceCellCentres() - Cf()))
        );

        w = nfc/(mag(n & (Cf() - Cn())) + nfc + SMALL);
    }
    else
    {
        // Slave side. Interpolate the master side weights, scale them for
        // partially covered faces and set weights for fully uncovered faces if
        // the bridge overlap is switched on. VV, 15/Feb/2018.

        // Pick up weights from the master side
        scalarField masterWeights(shadow().size());
        shadow().makeWeights(masterWeights);

        // Interpolate master weights to this side
        w = interpolate(masterWeights);

        if (bridgeOverlap())
        {
            // Weights for fully uncovered faces
            const scalarField uncoveredWeights(w.size(), 0.5);

            // Set weights for uncovered faces
            setUncoveredFaces(uncoveredWeights, w);

            // Scale partially overlapping faces
            scalePartialFaces(w);
        }

        // Finally construct these weights as 1 - master weights
        w = 1.0 - w;
    }
}


// Make patch face - neighbour cell distances
void Foam::ggiFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (ggiPolyPatch_.owner())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.

        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        const vectorField d(delta());

        dc = 1.0/max(nf() & d, 0.05*mag(d));

        // Note: no need to bridge the overlap since delta already takes it into
        // account. VV, 18/Oct/2017.
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltaCoeffs for fully uncovered faces if the
        // bridge overlap is switched on. VV, 15/Feb/2018.

        scalarField masterDeltas(shadow().size());
        shadow().makeDeltaCoeffs(masterDeltas);
        dc = interpolate(masterDeltas);

        if (bridgeOverlap())
        {
            // Delta coeffs for fully uncovered faces obtained from deltas on
            // this side
            const vectorField d(delta());
            const scalarField uncoveredDeltaCoeffs
            (
                1.0/max(nf() & d, 0.05*mag(d))
            );

            // Set delta coeffs for uncovered faces
            setUncoveredFaces(uncoveredDeltaCoeffs, dc);

            // Scale partially overlapping faces
            scalePartialFaces(dc);
        }
    }
}


// Make patch face non-orthogonality correction vectors
void Foam::ggiFvPatch::makeCorrVecs(vectorField& cv) const
{
    // Non-orthogonality correction on a ggi interface
    // MB, 7/April/2009

    // No non-orthogonal correction if the bridge overlap is switched on to
    // ensure conservative interpolation for partially overlapping faces
    if (bridgeOverlap())
    {
        cv = vector::zero;
    }
    else
    {
        // Calculate correction vectors on coupled patches
        const scalarField& patchDeltaCoeffs = deltaCoeffs();

        const vectorField patchDeltas(delta());
        const vectorField n(nf());

        // If non-orthogonality is over 90 deg, kill correction vector
        // HJ, 6/Jan/2011
        cv = pos(patchDeltas & n)*(n - patchDeltas*patchDeltaCoeffs);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::ggiFvPatch::delta() const
{
    if (ggiPolyPatch_.owner())
    {
        // Master side. Note: scaling partially covered faces and setting deltas
        // to fully uncovered faces correctly taken into account in
        // reconFaceCellCentres function. VV, 15/Feb/2018.

        tmp<vectorField> tdelta =
            ggiPolyPatch_.reconFaceCellCentres() - Cn();

        return tdelta;
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltas for fully uncovered faces if the bridge
        // overlap is switched on. VV, 15/Feb/2018.

        tmp<vectorField> tdelta = interpolate
        (
            shadow().Cn() - ggiPolyPatch_.shadow().reconFaceCellCentres()
        );

        if (bridgeOverlap())
        {
            // Deltas for fully uncovered faces
            const vectorField uncoveredDeltas(2.0*fvPatch::delta());

            // Set deltas for fully uncovered faces
            setUncoveredFaces(uncoveredDeltas, tdelta.ref());

            // Scale for partially covered faces
            scalePartialFaces(tdelta.ref());
        }

        return tdelta;
    }
}


const Foam::ggiFvPatch& Foam::ggiFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[ggiPolyPatch_.shadowIndex()];

    return refCast<const ggiFvPatch>(p);
}


bool Foam::ggiFvPatch::owner() const
{
    return ggiPolyPatch_.owner();
}


bool Foam::ggiFvPatch::fineLevel() const
{
    return true;
}


Foam::label Foam::ggiFvPatch::shadowIndex() const
{
    return ggiPolyPatch_.shadowIndex();
}


const Foam::ggiLduInterface& Foam::ggiFvPatch::shadowInterface() const
{
    const fvPatch& p = this->boundaryMesh()[ggiPolyPatch_.shadowIndex()];

    return refCast<const ggiLduInterface>(p);
}


Foam::label Foam::ggiFvPatch::interfaceSize() const
{
    return ggiPolyPatch_.size();
}


Foam::label Foam::ggiFvPatch::zoneSize() const
{
    return ggiPolyPatch_.zone().size();
}


const Foam::labelList& Foam::ggiFvPatch::zoneAddressing() const
{
    return ggiPolyPatch_.zoneAddressing();
}


const Foam::labelListList& Foam::ggiFvPatch::ggiAddressing() const
{
    if (ggiPolyPatch_.owner())
    {
        return ggiPolyPatch_.patchToPatch().masterAddr();
    }
    else
    {
        return ggiPolyPatch_.patchToPatch().slaveAddr();
    }
}


bool Foam::ggiFvPatch::localParallel() const
{
    return ggiPolyPatch_.localParallel();
}


const Foam::mapDistribute& Foam::ggiFvPatch::map() const
{
    return ggiPolyPatch_.map();
}


const Foam::scalarListList& Foam::ggiFvPatch::ggiWeights() const
{
    if (ggiPolyPatch_.owner())
    {
        return ggiPolyPatch_.patchToPatch().masterWeights();
    }
    else
    {
        return ggiPolyPatch_.patchToPatch().slaveWeights();
    }
}


void Foam::ggiFvPatch::expandAddrToZone(labelField& lf) const
{
    lf = ggiPolyPatch_.fastExpand(lf);
}


void Foam::ggiFvPatch::expandCrMatrixToZone(crMatrix& patchP) const
{
    if (!localParallel())
    {
        // Split the crMatrix into rows and expand it
        const crAddressing& patchCrAddr = patchP.crAddr();
        const labelList& patchRowStart = patchCrAddr.rowStart();
        const labelList& patchCol = patchCrAddr.column();
        const scalarField& patchCoeff = patchP.coeffs();

        List<labelField> cols(patchCrAddr.nRows());
        List<scalarField> coeffs(patchCrAddr.nRows());

        for (label faceI = 0; faceI < patchCrAddr.nRows(); faceI++)
        {
            // Unpack row
            const label rowStart = patchRowStart[faceI];
            const label rowLength = patchRowStart[faceI + 1] - rowStart;

            cols[faceI].setSize(rowLength);
            labelField& curCols = cols[faceI];

            coeffs[faceI].setSize(rowLength);
            scalarField& curCoeffs = coeffs[faceI];

            for (label coeffI = 0; coeffI < rowLength; coeffI++)
            {
                curCols[coeffI] = patchCol[rowStart + coeffI];
                curCoeffs[coeffI] = patchCoeff[rowStart + coeffI];
            }
        }

        // Expand to zone size
        List<labelField> zoneColsFF = ggiPolyPatch_.fastExpand(cols);
        List<scalarField> zoneCoeffsFF = ggiPolyPatch_.fastExpand(coeffs);

        scalar nZoneEntries = 0;

        forAll (zoneColsFF, zfI)
        {
            nZoneEntries += zoneColsFF[zfI].size();
        }

        // Reconstruct matrix
        labelList zoneRowStart(zoneSize() + 1);
        labelList zoneCols(nZoneEntries);
        scalarField zoneCoeffs(nZoneEntries);

        zoneRowStart[0] = 0;
        // Reset nZoneEntries for use as a counter
        nZoneEntries = 0;

        forAll(zoneColsFF, zfI)
        {
            const labelField& curCols = zoneColsFF[zfI];
            const scalarField& corCoeffs = zoneCoeffsFF[zfI];

            zoneRowStart[zfI + 1] = zoneRowStart[zfI] + curCols.size();

            forAll (curCols, coeffI)
            {
                zoneCols[nZoneEntries] = curCols[coeffI];
                zoneCoeffs[nZoneEntries] = corCoeffs[coeffI];
                nZoneEntries++;
            }
        }
        patchP = crMatrix
        (
            zoneSize(),
            patchCrAddr.nCols(),
            zoneRowStart,
            zoneCols
        );

        // Set coeffs
        patchP.coeffs() = zoneCoeffs;
    }
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::ggiFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::transfer
(
    const Pstream::commsTypes,
    const labelUList& interfaceData
) const
{
    return this->shadow().labelTransferBuffer();
}


void Foam::ggiFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    labelUList& iF
) const
{
    // Label transfer is local without global reduction
    labelTransferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    return shadow().labelTransferBuffer();
}


void Foam::ggiFvPatch::initProlongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    // crMatrix transfer is local without global reduction
    crMatrixTransferBuffer_ = filteredP;
}


Foam::autoPtr<Foam::crMatrix> Foam::ggiFvPatch::prolongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    autoPtr<crMatrix> tnbrP(new crMatrix(shadow().crMatrixTransferBuffer()));

    return tnbrP;
}


// ************************************************************************* //
