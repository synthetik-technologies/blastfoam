/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "fvMeshPrimitiveLduAddressing.H"
#include "lduPrimitiveMesh.H"
#include "processorLduInterface.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::fvMeshPrimitiveLduAddressing::fvMeshPrimitiveLduAddressing
(
    const fvMesh& mesh
)
:
    lduAddressing(mesh.nCells()),
    lowerAddr_
    (
        labelList::subList
        (
            mesh.faceOwner(),
            mesh.nInternalFaces()
        )
    ),
    upperAddr_(mesh.faceNeighbour()),
    patchAddr_(mesh.boundary().size()),
    patchSchedule_(mesh.globalData().patchSchedule())
{
    forAll(mesh.boundary(), patchI)
    {
        patchAddr_[patchI] = &mesh.boundary()[patchI].faceCells();
    }
}


Foam::fvMeshPrimitiveLduAddressing::fvMeshPrimitiveLduAddressing
(
    const label nCells,
    labelList&& lowerAddr,
    labelList&& upperAddr,
    const List<const labelUList*>& patchAddr,
    const lduSchedule& ps
)
:
    lduAddressing(nCells),
    lowerAddr_(std::move(lowerAddr)),
    upperAddr_(std::move(upperAddr)),
    patchAddr_(patchAddr),
    patchSchedule_(ps)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fvMeshPrimitiveLduAddressing::triIndex
(
    const lduAddressing& addr,
    const label a,
    const label b
)
{
    label own = min(a, b);

    label nbr = max(a, b);

    label startLabel = addr.ownerStartAddr()[own];

    label endLabel = addr.ownerStartAddr()[own + 1];

    const labelUList& neighbour = addr.upperAddr();

    for (label i = startLabel; i < endLabel; i++)
    {
        if (neighbour[i] == nbr)
        {
            return i;
        }
    }
    return -1;
}


Foam::labelList Foam::fvMeshPrimitiveLduAddressing::addAddressing
(
    const lduAddressing& addr,
    const labelListList& nbrCells,
    label& nExtraFaces,
    labelList& newLowerAddr,
    labelList& newUpperAddr,
    labelListList& nbrCellFaces,
    const globalIndex& globalNumbering,
    const labelList& globalCellIDs,
    labelListList& localFaceCells,
    labelListList& remoteFaceCells
)
{
    label nCells = addr.size();
    label nFaces = addr.upperAddr().size();
    labelList nProcFaces(Pstream::nProcs(), Zero);

    // Count additional faces
    nExtraFaces = 0;
    forAll(nbrCells, cellI)
    {
        const labelList& nbrs = nbrCells[cellI];
        forAll(nbrs, nbrI)
        {
            if (nbrs[nbrI] < nCells)
            {
                // Local cell
                if (triIndex(addr, cellI, nbrs[nbrI]) == -1)
                {
                    nExtraFaces++;
                }
            }
            else
            {
                label globalNbr = globalCellIDs[nbrs[nbrI]];
                label procI = globalNumbering.whichProcID(globalNbr);
                nProcFaces[procI]++;
            }
        }
    }

    // Create space for extra addressing
    newLowerAddr.setSize(nFaces + nExtraFaces);
    newUpperAddr.setSize(nFaces + nExtraFaces);

    // Copy existing addressing
    SubList<label>(newLowerAddr, nFaces) = addr.lowerAddr();
    SubList<label>(newUpperAddr, nFaces) = addr.upperAddr();


    // Per processor its local cells we want
    localFaceCells.setSize(Pstream::nProcs());
    remoteFaceCells.setSize(Pstream::nProcs());
    forAll(nProcFaces, procI)
    {
        localFaceCells[procI].setSize(nProcFaces[procI]);
        remoteFaceCells[procI].setSize(nProcFaces[procI]);
    }
    nProcFaces = 0;

    nbrCellFaces.setSize(nbrCells.size());
    forAll(nbrCells, cellI)
    {
        const labelList& nbrs = nbrCells[cellI];
        labelList& faces = nbrCellFaces[cellI];
        faces.setSize(nbrs.size());

        forAll(nbrs, nbrI)
        {
            label nbrCellI = nbrs[nbrI];

            if (nbrCellI < nCells)
            {
                // Find neighbour cell in owner list
                label faceI = triIndex(addr, cellI, nbrCellI);
                if (faceI == -1)
                {
                    faceI = nFaces++;
                    newLowerAddr[faceI] = min(cellI, nbrCellI);
                    newUpperAddr[faceI] = max(cellI, nbrCellI);
                }
                faces[nbrI] = faceI;
            }
            else
            {
                // Remote cell
                faces[nbrI] = -1;

                label globalNbr = globalCellIDs[nbrCellI];
                label procI = globalNumbering.whichProcID(globalNbr);
                label remoteCellI = globalNumbering.toLocal(procI, globalNbr);

                label procFaceI = nProcFaces[procI]++;
                localFaceCells[procI][procFaceI] = cellI;
                remoteFaceCells[procI][procFaceI] = remoteCellI;
            }
        }
    }

    // Reorder upper-triangular
    labelList oldToNew
    (
        upperTriOrder
        (
            addr.size(),
            newLowerAddr,
            newUpperAddr
        )
    );

    // Shuffle face-to-cell addressing
    inplaceReorder(oldToNew, newLowerAddr);
    inplaceReorder(oldToNew, newUpperAddr);
    // Update cell-to-face addressing
    forAll(nbrCellFaces, cellI)
    {
        inplaceRenumber(oldToNew, nbrCellFaces[cellI]);
    }

    //if (debug)
    //{
    //    for
    //    (
    //        label i = addr.upperAddr().size();
    //        i < oldToNew.size();
    //        i++
    //    )
    //    {
    //        label faceI = oldToNew[i];
    //        Pout<< "new face:" << faceI << endl
    //            << "    owner:" << newLowerAddr[faceI]
    //            << "    neighbour:" << newUpperAddr[faceI]
    //            << endl;
    //    }
    //
    //    forAll(nbrCellFaces, cellI)
    //    {
    //        const labelList& nbrs = nbrCells[cellI];
    //        const labelList& nbrFaces = nbrCellFaces[cellI];
    //        if (nbrs.size())
    //        {
    //            Pout<< "cell:" << cellI << " has additional neighbours:"
    //                << endl;
    //            forAll(nbrs, i)
    //            {
    //                label faceI = nbrFaces[i];
    //                Pout<< "    nbr:" << nbrs[i]
    //                    << " through face:" << faceI
    //                    << " with own:" << newLowerAddr[faceI]
    //                    << " with nei:" << newUpperAddr[faceI]
    //                    << endl;
    //            }
    //        }
    //    }
    //}

    return oldToNew;
}


Foam::labelList Foam::fvMeshPrimitiveLduAddressing::upperTriOrder
(
    const label nCells,
    const labelUList& lower,
    const labelUList& upper
)
{
    labelList nNbrs(nCells, 0);

    // Count number of upper neighbours
    forAll(lower, facei)
    {
        if (upper[facei] < lower[facei])
        {
            FatalErrorInFunction
                << "Problem at face:" << facei
                << " lower:" << lower[facei]
                << " upper:" << upper[facei]
                << exit(FatalError);
        }
        nNbrs[lower[facei]]++;
    }

    // Construct cell-upper cell addressing
    labelList offsets(nCells+1);
    offsets[0] = 0;
    forAll(nNbrs, celli)
    {
        offsets[celli+1] = offsets[celli]+nNbrs[celli];
    }

    nNbrs = offsets;

    labelList cellToFaces(offsets.last());
    forAll(upper, facei)
    {
        label celli = lower[facei];
        cellToFaces[nNbrs[celli]++] = facei;
    }

    // Sort

    labelList oldToNew(lower.size());

    labelList order;
    labelList nbr;

    label newFacei = 0;

    for (label celli = 0; celli < nCells; celli++)
    {
        label startOfCell = offsets[celli];
        label nNbr = offsets[celli+1] - startOfCell;

        nbr.setSize(nNbr);
        order.setSize(nNbr);
        forAll(order, i)
        {
            nbr[i] = upper[cellToFaces[offsets[celli]+i]];
        }
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            label index = order[i];
            oldToNew[cellToFaces[startOfCell + index]] = newFacei++;
        }
    }

    return oldToNew;
}

// ************************************************************************* //
