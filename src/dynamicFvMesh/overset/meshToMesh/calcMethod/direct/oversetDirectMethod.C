/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "oversetDirectMethod.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetDirectMethod, 0);
    addToRunTimeSelectionTable(oversetMeshToMeshMethod, oversetDirectMethod, meshComponents);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::oversetDirectMethod::intersect
(
    const label srcCelli,
    const label tgtCelli
) const
{
    return tgt_.pointInCell
    (
        src_.cellCentres()[srcCelli],
        tgtCelli,
        polyMesh::FACE_PLANES
    );
}


bool Foam::oversetDirectMethod::findInitialSeeds
(
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const label startSeedI,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const cellList& srcCells = src_.cells();
    const faceList& srcFaces = src_.faces();
    const pointField& srcPts = src_.points();

    for (label i = startSeedI; i < srcCellIDs.size(); i++)
    {
        label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const point srcCtr(srcCells[srcI].centre(srcPts, srcFaces));
            label tgtI = tgt_.cellTree().findInside(srcCtr);

            if (tgtI != -1 && intersect(srcI, tgtI))
            {
                srcSeedI = srcI;
                tgtSeedI = tgtI;

                return true;
            }
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


void Foam::oversetDirectMethod::calculateAddressing
(
    labelListList& srcToTgtCellAddr,
    scalarListList& srcToTgtCellWght,
    labelListList& tgtToSrcCellAddr,
    scalarListList& tgtToSrcCellWght,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs, // not used
    boolList& mapFlag,
    label& startSeedI
)
{
    // store a list of src cells already mapped
    labelList srcTgtSeed(src_.nCells(), -1);

    List<DynamicList<label>> srcToTgt(src_.nCells());
    List<DynamicList<label>> tgtToSrc(tgt_.nCells());

    DynamicList<label> srcSeeds(10);

    const scalarField& srcVc = src_.cellVolumes();
    const scalarField& tgtVc = tgt_.cellVolumes();

    label srcCelli = srcSeedI;
    label tgtCelli = tgtSeedI;

    do
    {
        // store src/tgt cell pair
        srcToTgt[srcCelli].append(tgtCelli);
        tgtToSrc[tgtCelli].append(srcCelli);

        // mark source cell srcSeedI as matched
        mapFlag[srcCelli] = false;

        // accumulate intersection volume
        V_ += srcVc[srcCelli];

        // find new source seed cell
        appendToDirectSeeds
        (
            mapFlag,
            srcTgtSeed,
            srcSeeds,
            srcCelli,
            tgtCelli
        );
    }
    while (srcCelli >= 0);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellWght[i] = scalarList(srcToTgt[i].size(), srcVc[i]);
        srcToTgtCellAddr[i].transfer(srcToTgt[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellWght[i] = scalarList(tgtToSrc[i].size(), tgtVc[i]);
        tgtToSrcCellAddr[i].transfer(tgtToSrc[i]);
    }
}


void Foam::oversetDirectMethod::appendToDirectSeeds
(
    boolList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const labelList& srcNbr = src_.cellCells()[srcSeedI];
    const labelList& tgtNbr = tgt_.cellCells()[tgtSeedI];

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if (mapFlag[srcI] && (srcTgtSeed[srcI] == -1))
        {
            // source cell srcI not yet mapped

            // identify if target cell exists for source cell srcI
            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];

                if (intersect(srcI, tgtI))
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            if (!found)
            {
                // no match available for source cell srcI
                mapFlag[srcI] = false;
            }
        }
    }

    if (srcSeeds.size())
    {
        srcSeedI = srcSeeds.remove();
        tgtSeedI = srcTgtSeed[srcSeedI];
    }
    else
    {
        srcSeedI = -1;
        tgtSeedI = -1;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetDirectMethod::oversetDirectMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    oversetMeshToMeshMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetDirectMethod::~oversetDirectMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oversetDirectMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    List<List<point>>& srcToTgtVec,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght,
    List<List<point>>& tgtToSrcVec
)
{
    bool ok = initialise
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    if (!ok)
    {
        return;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells());

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src_.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        calculateAddressing
        (
            srcToTgtAddr,
            srcToTgtWght,
            tgtToSrcAddr,
            tgtToSrcWght,
            srcSeedI,
            tgtSeedI,
            srcCellIDs,
            mapFlag,
            startSeedI
        );
    }
    else
    {
        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }
}


// ************************************************************************* //
