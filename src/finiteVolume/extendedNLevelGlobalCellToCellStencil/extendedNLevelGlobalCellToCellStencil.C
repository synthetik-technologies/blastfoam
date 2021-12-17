/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "extendedNLevelGlobalCellToCellStencil.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class StencilType>
void Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::unionEqOp::
operator()
(
    labelList& x,
    const labelList& y
) const
{
    if (y.size())
    {
        if (x.empty())
        {
            x = y;
        }
        else
        {
            labelHashSet set(x);
            forAll(y, i)
            {
                set.insert(y[i]);
            }
            x = set.toc();
        }
    }
}


template<class StencilType>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::addCellNeighbours
(
    const Map<labelList>& cellCells,
    const label celli,
    Map<label>& levels,
    DynamicList<label>& neighbours,
    const label level
) const
{
    if (level > nLevels_)
    {
        return;
    }

    // Check if this cell has been added yet
    Map<label>::iterator iter = levels.find(celli);

    if(iter == levels.end())
    {
        // Add the cell because it does not exist yet
        levels.insert(celli, level);
    }
    else if (iter() <= level)
    {
        // The current cell has been added at a lower level
        // so its neighbors have also been added
        return;
    }
    else
    {
        // This cell is being added at the lowest level
        levels[celli] = level;
    }

    if (level == nLevels_)
    {
        return;
    }

    //- Append all neighbours to the list (only unique)
    unionEqOp()(neighbours, cellCells[celli]);


    // Extend based on neighbor neighbors
    // Starts at 1 since this cell is 0
    const labelList& cc(cellCells[celli]);
    for (label i = 1; i < cc.size(); i++)
    {
        addCellNeighbours
        (
            cellCells,
            cc[i],
            levels,
            neighbours,
            level + 1
        );
    }
}


template<class StencilType>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::makeGlobal
(
    const globalIndex& oldGlobalIndex,
    labelList& indices
) const
{
    label ni = 0;
    forAll(indices, i)
    {
        label proci = oldGlobalIndex.whichProcID(indices[i]);
        label idx = oldGlobalIndex.toLocal(proci, indices[i]);
        if (idx < gIndexPtr_().localSize(proci))
        {
            indices[ni++] = gIndexPtr_().toGlobal(proci, idx);
        }
    }
    indices.resize(ni);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class StencilType>
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::
extendedNLevelGlobalCellToCellStencil
(
    const polyMesh& mesh,
    const label nLevels
)
:
    MeshObject
    <
        polyMesh,
        PatchMeshObject,
        extendedNLevelGlobalCellToCellStencil<StencilType>
    >(mesh),
    mesh_(mesh),
    nLevels_(nLevels),
    gIndexPtr_(),
    stencilMap_(),
    cellCells_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class StencilType>
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::
~extendedNLevelGlobalCellToCellStencil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class StencilType>
bool Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::movePoints()
{
    mapPtr_.clear();
    gIndexPtr_.clear();
    return true;
}


template<class StencilType>
void Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::updateMesh
(
    const mapPolyMesh& mpm
)
{
    mapPtr_.clear();
    gIndexPtr_.clear();
}


template<class StencilType>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    // Assuming balancing
    mapPtr_.clear();
    gIndexPtr_.clear();
}


template<class StencilType>
void Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::addPatch
(
    const label patchi
)
{
    // Assuming balancing
    mapPtr_.clear();
    gIndexPtr_.clear();
}


template<class StencilType>
Foam::autoPtr<Foam::mapDistribute> Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::buildMap
(
    const List<label>& toProc
)
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(toProc, i)
    {
        label proci = toProc[i];

        nSend[proci]++;
    }


    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());

    forAll(nSend, proci)
    {
        sendMap[proci].setSize(nSend[proci]);

        nSend[proci] = 0;
    }

    // 3. Fill sendMap
    forAll(toProc, i)
    {
        label proci = toProc[i];

        sendMap[proci][nSend[proci]++] = i;
    }

    // 4. Send over how many I need to receive
    labelList recvSizes;
    Pstream::exchangeSizes(sendMap, recvSizes);


    // Determine receive map
    // ~~~~~~~~~~~~~~~~~~~~~

    labelListList constructMap(Pstream::nProcs());

    // Local transfers first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label constructSize = constructMap[Pstream::myProcNo()].size();

    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            label nRecv = recvSizes[proci];

            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[proci][i] = constructSize++;
            }
        }
    }

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            constructSize,
            move(sendMap),
            move(constructMap)
        )
    );
}


template<class StencilType>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::updateStencil() const
{
    if (!gIndexPtr_.valid())
    {
        gIndexPtr_.set(new globalIndex(mesh_.nCells()));
    }
    cellCells_.resize(mesh_.nCells());
    stencilMap_.clear();

    Map<labelList> singleLevelMap;

    //- Create global map (1 level deep)
    {
        StencilType ctcStencil(mesh_);
        List<List<label>> stencil(ctcStencil);
        forAll(stencil, i)
        {
            makeGlobal(ctcStencil.globalNumbering(), stencil[i]);
            singleLevelMap.insert(gIndexPtr_->toGlobal(i), stencil[i]);
        }
    }

    //- Add additional cells across processors that are not owned
    //  by the current processor
    if (Pstream::parRun())
    {
        labelHashSet missingCells;
        forAllConstIter
        (
            Map<labelList>,
            singleLevelMap,
            iter
        )
        {
            forAll(iter(), j)
            {
                label cellj = iter()[j];
                if (!gIndexPtr_->isLocal(cellj))
                {
                    missingCells.insert(cellj);
                }
            }
        }

        while (returnReduce(missingCells.size(), sumOp<label>()))
        {
            // Create the requests for new cell neighbours
            DynamicList<label> requests(Pstream::nProcs());
            labelList requestedCells(missingCells.size());
            label idx = 0;
            forAllConstIter
            (
                labelHashSet,
                missingCells,
                iter
            )
            {
                label celli = iter.key();
                requests.append(gIndexPtr_->whichProcID(celli));
                requestedCells[idx++] = celli;
            }

            // Needed for reverseDistribute
            label constructSize = requests.size();

            // make the map
            autoPtr<mapDistribute> map(buildMap(requests));

            // Send requests
            map().distribute(requestedCells);

            // Add neighbours to be sent using the ordering provided
            labelListList newCells(requestedCells.size());
            forAll(requestedCells, i)
            {
                newCells[i] = singleLevelMap[requestedCells[i]];
            }

            // Send the data back
            map().reverseDistribute(constructSize, newCells);

            // Insert the new neighbours to the single level map
            forAll(newCells, i)
            {
                singleLevelMap.insert
                (
                    newCells[i][0],
                    newCells[i]
                );
            }

            // Update missing cells that were added with the new cells
            missingCells.clear();
            forAll(newCells, i)
            {
                forAll(newCells[i], j)
                {
                    if (!singleLevelMap.found(newCells[i][j]))
                    {
                        missingCells.insert(newCells[i][j]);
                    }
                }
            }
        }
    }

    // Build the expanding cellCell map
    forAll(cellCells_, celli)
    {
        Map<label> levels;
        DynamicList<label> neighbors;

        const label gCelli = gIndexPtr_->toGlobal(celli);
        addCellNeighbours
        (
            singleLevelMap,
            gCelli,
            levels,
            neighbors
        );
        stencilMap_.insert
        (
            gCelli,
            cellStencil
            (
                gCelli,
                move(neighbors),
                levels,
                mesh_.cellCentres()[celli]
            )
        );
        cellCells_[celli] = stencilMap_[gCelli];
    }

    // Collect cell stencils not on this processor
    if (Pstream::parRun())
    {
        Map<label> idxMap;
        List<labelList> sendMap(Pstream::nProcs());
        List<cellStencil> newCells;
        label idx = 0;

        // Check if local stencil
        forAll(cellCells_, i)
        {
            label celli = gIndexPtr_->toGlobal(i);
            boolList added(Pstream::nProcs(), false);
            bool oneAdded = false;
            const cellStencil& lStencil(stencilMap_[celli]);

            forAll(lStencil, j)
            {
                label cellj = lStencil[j];
                label proci = gIndexPtr_->whichProcID(cellj);
                if (proci != Pstream::myProcNo() && !added[proci])
                {
                    sendMap[proci].append(idx);
                    added[proci] = true;
                    oneAdded = true;
                }
            }
            if (oneAdded)
            {
                newCells.append(lStencil);
                idxMap.insert(celli, idx++);
            }
        }

        labelList recvSizes;
        Pstream::exchangeSizes(sendMap, recvSizes);

        // Determine receive map
        // ~~~~~~~~~~~~~~~~~~~~~
        labelListList constructMap(Pstream::nProcs());

        // Local transfers first
        constructMap[Pstream::myProcNo()] = identity
        (
            sendMap[Pstream::myProcNo()].size()
        );

        label constructSize =
            constructMap[Pstream::myProcNo()].size();

        forAll(constructMap, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                label nRecv = recvSizes[proci];

                constructMap[proci].setSize(nRecv);

                for (label i = 0; i < nRecv; i++)
                {
                    constructMap[proci][i] = constructSize++;
                }
            }
        }

        mapDistribute map
        (
            constructSize,
            move(sendMap),
            move(constructMap)
        );

        map.distribute(newCells);

        forAll(newCells, i)
        {
            stencilMap_.insert(newCells[i].owner(), newCells[i]);
        }
    }

    //- Clean up map
    //  remove stencils that do not have any local indicies,
    //  and convert to local
    forAllIter
    (
        Map<cellStencil>,
        stencilMap_,
        iter
    )
    {
        iter().updateLocalStencil(gIndexPtr_());
        if (!iter().localStencil().size())
        {
            stencilMap_.erase(iter);
        }
    }


    List<Map<label>> compactMap(Pstream::nProcs());
    mapPtr_.reset
    (
        new mapDistribute
        (
            gIndexPtr_(),
            cellCells_,
            compactMap
        )
    );
}


// ************************************************************************* //
