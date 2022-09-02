/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "waveMethod.H"
#include "meshToMeshData.H"
#include "FaceCellWave.H"
#include "addToRunTimeSelectionTable.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveMethod, 0);
    addToRunTimeSelectionTable(oversetMeshToMeshMethod, waveMethod, meshComponents);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::waveMethod::calculate
(
    const polyMesh& src,
    const polyMesh& tgt,
    labelList& srcToTgtAddr
)
{
    // If parallel running a local domain might have zero cells thus never
    // constructing the face-diagonal decomposition which uses parallel
    // transfers.
    (void)tgt.tetBasePtIs();

    // The actual matching is only w.r.t local cells so cannot be run in
    // parallel.
    const bool oldParRun = Pstream::parRun();
    Pstream::parRun() = false;

    label nSeeds = 0;

    if (tgt.nCells() == 0)
    {
        srcToTgtAddr.setSize(src.nCells());
        srcToTgtAddr = -1;
    }
    else
    {
        const treeBoundBox& tgtBb = tgt.cellTree().bb();

        DynamicList<label> changedFaces(src.nFaces()/100 + 100);
        DynamicList<meshToMeshData> changedFacesInfo(changedFaces.size());

        List<meshToMeshData> cellData(src.nCells());
        List<meshToMeshData> faceData(src.nFaces());

        meshToMeshData::trackData td(tgt);

        label startCelli = 0;

        while (true)
        {
            changedFaces.clear();
            changedFacesInfo.clear();

            // Search for starting seed
            for (; startCelli < src.nCells(); startCelli++)
            {
                if (!cellData[startCelli].valid(td))
                {
                    nSeeds++;
                    const point& cc = src.cellCentres()[startCelli];

                    if (!tgtBb.contains(cc))
                    {
                        // Point outside local bb of tgt mesh. No need to
                        // search. Register as no correspondence
                        cellData[startCelli] = meshToMeshData(-1);
                    }
                    else
                    {
                        label tgtCelli = tgt.findCell(cc, polyMesh::CELL_TETS);
                        if (tgtCelli != -1)
                        {
                            // Insert any face of cell
                            label facei = src.cells()[startCelli][0];
                            changedFaces.append(facei);
                            changedFacesInfo.append(meshToMeshData(tgtCelli));
                            break;
                        }
                        else
                        {
                            // Register as no correspondence
                            cellData[startCelli] = meshToMeshData(-1);
                        }
                    }
                }
            }

            if (returnReduce(changedFaces.empty(), andOp<bool>()))
            {
                break;
            }

            FaceCellWave<meshToMeshData, meshToMeshData::trackData> calc
            (
                src,
                changedFaces,
                changedFacesInfo,
                faceData,
                cellData,
                src.globalData().nTotalCells(),   // max iterations
                td
            );
        }

        // Copy into srcToTgt
        srcToTgtAddr.setSize(src.nCells());

        forAll(cellData, celli)
        {
            srcToTgtAddr[celli] = cellData[celli].tgtCell();
        }
    }

    Pstream::parRun() = oldParRun;  // Restore parallel state

    if (debug)
    {
        Pout<< "nSeeds:" << returnReduce(nSeeds, sumOp<label>())
            << " out of nCells:" << returnReduce(src.nCells(), sumOp<label>())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMethod::waveMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    oversetMeshToMeshMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveMethod::~waveMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    List<List<point>>& srcToTgtVec,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght,
    List<List<point>>& tgtToSrcVec
)
{
    {
        labelList srcToTgt(src_.nCells());
        calculate(src_, tgt_, srcToTgt);
        srcToTgtAddr.setSize(srcToTgt.size());
        srcToTgtWght.setSize(srcToTgt.size());
        forAll(srcToTgtAddr, celli)
        {
            srcToTgtAddr[celli].setSize(1);
            srcToTgtAddr[celli][0] = srcToTgt[celli];
            srcToTgtWght[celli].setSize(1);
            srcToTgtWght[celli][0] = src_.cellVolumes()[celli];
        }
    }

    {
        labelList tgtToSrc(tgt_.nCells());
        calculate(tgt_, src_, tgtToSrc);
        tgtToSrcAddr.setSize(tgtToSrc.size());
        tgtToSrcWght.setSize(tgtToSrc.size());
        forAll(tgtToSrcAddr, celli)
        {
            tgtToSrcAddr[celli].setSize(1);
            tgtToSrcAddr[celli][0] = tgtToSrc[celli];
            tgtToSrcWght[celli].setSize(1);
            tgtToSrcWght[celli][0] = tgt_.cellVolumes()[celli];
        }
    }
}


// ************************************************************************* //
