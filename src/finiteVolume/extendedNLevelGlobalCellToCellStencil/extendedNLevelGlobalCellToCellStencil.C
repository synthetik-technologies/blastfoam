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
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::addCellNeighbours
(
    const Map<labelList>& cellCells,
    const label celli,
    const label level,
    const label maxLevel,
    labelHashSet& neighbours,
    labelHashSet& added,
    DynamicList<label>& nbrs
) const
{
	Map<labelList>::const_iterator iterI = cellCells.find(celli);
	if (level > maxLevel || iterI == cellCells.cend())
	{
	    return;
	}
	const labelList& cc(iterI());

    // First pass, add neighours of neighbours
    forAll(cc, i)
    {
    	Map<labelList>::const_iterator iter = cellCells.find(cc[i]);
    	if (iter != cellCells.cend())
    	{
            forAll(iter(), j)
            {
                if (neighbours.insert(iter()[j]))
                {
                    nbrs.append(iter()[j]);
                }
            }
		}
    }

    if (level+1 == maxLevel)
    {
    	return;
    }

    // Extend based on neighbor neighbors
    forAll(cc, i)
    {
    	if (added.insert(cc[i]))
    	{
	        addCellNeighbours
	        (
	            cellCells,
	            cc[i],
	            level + 1,
	            maxLevel,
	            neighbours,
	            added,
                nbrs
	        );
        }
    }
}


template<class StencilType>
Foam::labelList
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::addCellNeighbours
(
    const Map<labelList>& cellCells,
    const label celli,
    labelHashSet& neighbours,
    DynamicList<label>& nbrs
) const
{
	Map<labelList>::const_iterator iterI = cellCells.find(celli);
	if (iterI == cellCells.cend())
	{
	    return {};
	}

	const labelList& cc(iterI());
	DynamicList<label> newCells(cc.size());

    // Add neighbours and mark new cells to be added
    forAll(cc, i)
    {
    	Map<labelList>::const_iterator iter = cellCells.find(cc[i]);
    	if (iter != cellCells.cend())
    	{
            forAll(iter(), j)
            {
                if (neighbours.insert(iter()[j]))
                {
                    nbrs.append(iter()[j]);
                }
            }
			newCells.append(cc[i]);
		}
    }

    return move(newCells);
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


template<class StencilType>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::createMap() const
{
	cellCells_.resize(this->mesh_.nCells());
	const globalIndex& gIndex = this->gIndex();
	forAllConstIter(Map<cellStencil>, cellCellMap(), iter)
	{
		if (gIndex.isLocal(iter.key()))
		{
			cellCells_[gIndex.toLocal(iter.key())] = iter();
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


template<class StencilType>
const Foam::mapDistribute&
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::map() const
{
	if (!mapPtr_.valid())
	{
		createMap();
	}
	return mapPtr_();
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
        UpdateableMeshObject,
        extendedNLevelGlobalCellToCellStencil<StencilType>
    >(mesh),
    mesh_(mesh),
    nLevels_(nLevels),
    nNbrs_(),
    gIndexPtr_(),
    stencilMap_(),
    cellCells_(),
    mapPtr_(),
    nonlocalCells_(),
    nonlocalOwners_()
{}


template<class StencilType>
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::
extendedNLevelGlobalCellToCellStencil
(
    const polyMesh& mesh,
    const labelList& nNbrs
)
:
    MeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        extendedNLevelGlobalCellToCellStencil<StencilType>
    >(mesh),
    mesh_(mesh),
    nLevels_(-1),
    nNbrs_(nNbrs),
    gIndexPtr_(),
    stencilMap_(),
    cellCells_(),
    mapPtr_(),
    nonlocalCells_(),
    nonlocalOwners_()
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
    return true;
}


template<class StencilType>
void Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::updateMesh
(
    const mapPolyMesh& mpm
)
{
	stencilMap_.clear();
    mapPtr_.clear();
    gIndexPtr_.clear();
    nonlocalCells_.clear();
    nonlocalOwners_.clear();
}


template<class StencilType>
Foam::autoPtr<Foam::mapDistribute>
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::buildMap
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
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::calcStencil() const
{
    if (!gIndexPtr_.valid())
    {
        gIndexPtr_.set(new globalIndex(mesh_.nCells()));
    }

    Map<labelList> singleLevelMap;

    //- Create global map (1 level deep)
    {
        StencilType ctcStencil(mesh_);
        List<List<label>> stencil(ctcStencil);
        forAll(stencil, i)
        {
            makeGlobal(ctcStencil.globalNumbering(), stencil[i]);
            singleLevelMap.insert
            (
            	gIndexPtr_->toGlobal(i),
            	stencil[i]
        	);
        }
    }

    //- Add additional cells across processors that are not owned
    //  by the current processor
    if (Pstream::parRun())
    {
        labelHashSet missingCells;
        label n = 0;
        forAllConstIter
        (
            Map<labelList>,
            singleLevelMap,
            iter
        )
        {
        	n += iter().size();
            forAll(iter(), j)
            {
                label cellj = iter()[j];
                if (!gIndexPtr_->isLocal(cellj))
                {
                    missingCells.insert(cellj);
                    const label proci = gIndexPtr_->whichProcID(cellj);
                    nonlocalCells_.insert(cellj, proci);
                }
            }
        }

        label leveli = 1;
        const label nTotal = gSum(nNbrs_);
        while (true)
        {
        	if (nNbrs_.size() == mesh_.nCells())
	        {
	        	if (returnReduce(n, sumOp<label>()) >= nTotal)
		        {
		        	break;
		        }
	        }
	        else if (nLevels_ > 0)
	        {
	        	if (leveli++ > nLevels_)
	        	{
	        		break;
	        	}
	        }

	        if (!missingCells.size())
	        {
	        	break;
	        }


            // Create the requests for new cell neighbours
            DynamicList<label> requests(missingCells.size());
            DynamicList<label> requestedCells(missingCells.size());
            forAllConstIter
            (
                labelHashSet,
                missingCells,
                iter
            )
            {
                label celli = iter.key();
                requests.append(gIndexPtr_->whichProcID(celli));
                requestedCells.append(celli);
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
            map().reverseDistribute(constructSize, requestedCells);

            // Insert the new neighbours to the single level map
            forAll(newCells, i)
            {
                singleLevelMap.insert
                (
                    requestedCells[i],
                    newCells[i]
                );
                n += newCells.size();
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
    forAll(this->mesh_.cells(), celli)
    {
        const label gCelli = gIndexPtr_->toGlobal(celli);
        labelHashSet neighbors(singleLevelMap[gCelli]);
        DynamicList<label> nbrs(16);
        labelHashSet added({gCelli});
        if (nNbrs_.size() == mesh_.nCells())
        {
        	labelList cc(singleLevelMap[gCelli]);
            nbrs.append(cc);
        	while (neighbors.size() < nNbrs_[celli] && cc.size())
        	{
        		labelHashSet newCells(cc.size());
        		forAll(cc, i)
        		{
        			if (added.insert(cc[i]))
        			{
	        			newCells.insert
	        			(
	        				addCellNeighbours
					        (
					            singleLevelMap,
					            cc[i],
					            neighbors,
                                nbrs
					        )
				        );
			        }
			    }
                newCells -= added;
			    cc = newCells.toc();
        	}
        }
        else if (nLevels_ > 0)
        {
            nbrs.append(singleLevelMap[gCelli]);
        	label leveli = 2;
        	addCellNeighbours
	        (
	            singleLevelMap,
	            gCelli,
	            leveli,
	            nLevels_,
	            neighbors,
	            added,
                nbrs
	        );
        }

        stencilMap_.insert
        (
            gCelli,
            cellStencil
            (
                gCelli,
                nbrs,
                mesh_.cellCentres()[celli]
            )
        );
    }

    update();
}


template<class StencilType>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::update() const
{
	cellCellMap();
	if (!Pstream::parRun())
	{
		return;
	}

	{
		// Create the requests for new cell neighbours
	    DynamicList<label> requests(stencilMap_.size());
	    DynamicList<cellStencil> sendCells(stencilMap_.size());
        boolList marked(Pstream::nProcs());
	    forAllConstIter(Map<cellStencil>, stencilMap_, iter)
	    {
            const labelList& cells = iter();
            marked = false;
            forAll(cells, cj)
            {
                const label procj = gIndexPtr_->whichProcID(cells[cj]);
                if (procj != Pstream::myProcNo() && !marked[procj])
                {
        	        requests.append(procj);
        	        sendCells.append(iter());
                    marked[procj] = true;
                }
            }
	    }

	    // make the map
	    autoPtr<mapDistribute> map(buildMap(requests));

	    // Send requests
	    map().distribute(sendCells);

	    // Add neighbours to be sent using the ordering provided
	    forAll(sendCells, i)
	    {
	        stencilMap_.insert(sendCells[i].owner(), move(sendCells[i]));
	    }
	}

	// //- Clean up map
    nonlocalCells_.clear();
    nonlocalOwners_.clear();
    forAllIter
    (
        Map<cellStencil>,
        stencilMap_,
        iter
    )
    {
        iter().update(gIndexPtr_());
        const label proci = gIndexPtr_->whichProcID(iter().owner());
        if (!iter().isLocal())
        {
            nonlocalOwners_.set(iter().owner(), proci);
        }

        const labelList& stencil = iter();
        forAll(stencil, i)
        {
        	const label procj = gIndexPtr_->whichProcID(stencil[i]);
        	if (procj != Pstream::myProcNo())
        	{
            	nonlocalCells_.set(stencil[i], procj);
        	}
        }
    }
}
// ************************************************************************* //
