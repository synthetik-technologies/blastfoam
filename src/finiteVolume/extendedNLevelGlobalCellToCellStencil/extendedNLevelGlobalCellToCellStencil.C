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
#include "extendedCentredCellToCellStencil.H"
#include "CECCellToCellStencil.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedNLevelGlobalCellToCellStencil, 0);
}

void Foam::cellStencil::updateLocalStencil
(
    const globalIndex& idx
)
{
    const labelList& stencil(*this);
    localStencil_ = stencil;
    label ni = 0;
    forAll(stencil, i)
    {
        if (idx.isLocal(stencil[i]))
        {
            localStencil_[ni++] = idx.toLocal(stencil[i]);
        }
    }
    localStencil_.resize(ni);
}

Foam::Ostream& Foam::operator<<(Ostream& os, const cellStencil& c)
{
    const labelList& lst(c);
    os << lst << " "
        << c.owner() << " "
        << c.centre() << " ";
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, cellStencil& c)
{
    labelList& lst(c);
    is >> lst;
    is >> c.owner();
    is >> c.centre();
    return is;
}

void Foam::extendedNLevelGlobalCellToCellStencil::unionEqOp::operator()
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

void Foam::extendedNLevelGlobalCellToCellStencil::addCellNeighbours
(
    const Map<labelList>& cellCells,
    const label celli,
    Map<label>& levels,
    labelList& neighbours,
    const label level
) const
{
    if (level > nLevels_ || !cellCells.found(celli))
    {
        return;
    }
    // Set the level
    if (levels.found(celli))
    {
        if (levels[celli] <= level)
        {
            return;
        }
        levels[celli] = level;
    }
    else
    {
        levels.insert(celli, level);
    }

    //- Append all neighbours to the list
    const labelList& cc(cellCells[celli]);
    unionEqOp()(neighbours, cc);

    // Extend based on neighbor neighbors
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

void Foam::extendedNLevelGlobalCellToCellStencil::removeGlobalFaces
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
        if (idx < gIndex_.localSize(proci))
        {
            indices[ni++] = gIndex_.toGlobal(proci, idx);
        }
    }
    indices.resize(ni);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedNLevelGlobalCellToCellStencil::extendedNLevelGlobalCellToCellStencil
(
    const polyMesh& mesh,
    const label nLevels
)
:
    NLevelStencil(mesh),
    mesh_(mesh),
    nLevels_(nLevels),
    gIndex_(mesh.nCells()),
    stencilMap_(),
    cellCells_(mesh.nCells())
{
    updateStencil();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedNLevelGlobalCellToCellStencil::~extendedNLevelGlobalCellToCellStencil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedNLevelGlobalCellToCellStencil::updateStencil() const
{
    gIndex_ = globalIndex(mesh_.nCells());
    cellCells_.resize(mesh_.nCells());
    stencilMap_.clear();

    Map<labelList> singleLevelMap;

    //- Create global map (1 level deep)
    {
        CECCellToCellStencil ctcStencil(mesh_);
        List<List<label>> stencil(ctcStencil);
        const globalIndex& cfgIndex = ctcStencil.globalNumbering();
        forAll(stencil, i)
        {
            removeGlobalFaces(cfgIndex, stencil[i]);
            singleLevelMap.insert(gIndex_.toGlobal(i), stencil[i]);
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
                if (!gIndex_.isLocal(cellj))
                {
                    missingCells.insert(cellj);
                }
            }
        }

        for (label leveli = 1; leveli < nLevels_; leveli++)
        {
            // Create the requests for new cell neighbours
            List<labelList> requests(Pstream::nProcs());
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
                requests[gIndex_.whichProcID(celli)].append(idx);
                requestedCells[idx++] = celli;

            }

            // Count number of requests
            labelList nSend(Pstream::nProcs());
            forAll(nSend, proci)
            {
                nSend[proci] = requests[proci].size();
            }

            labelList recvSizes;
            Pstream::exchangeSizes(requests, recvSizes);

            // Determine receive map
            // ~~~~~~~~~~~~~~~~~~~~~
            labelListList constructMap(Pstream::nProcs());

            // Local transfers first
            constructMap[Pstream::myProcNo()] = identity
            (
                requests[Pstream::myProcNo()].size()
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

            // Create map
            mapDistribute map
            (
                constructSize,
                move(requests),
                move(constructMap)
            );

            // Send requests
            map.distribute(requestedCells);

            // Add neighbours to be sent using the ordering provided
            labelListList newCells(requestedCells.size());
            forAll(requestedCells, i)
            {
                newCells[i] = singleLevelMap[requestedCells[i]];
            }

            // Send the data back
            map.reverseDistribute(constructSize, newCells);

            // Insert the new neighbours to the single level map
            forAll(newCells, i)
            {
                singleLevelMap.insert
                (
                    newCells[i][0],
                    newCells[i]
                );
            }
            if (leveli != nLevels_ - 1)
            {
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
    }

    // Build the expanding cellCell map
    Map<label> levels;
    forAll(cellCells_, celli)
    {
        levels.clear();

        labelList neighbors;
        addCellNeighbours
        (
            singleLevelMap,
            gIndex_.toGlobal(celli),
            levels,
            neighbors
        );

        const label gCelli = gIndex_.toGlobal(celli);
        stencilMap_.insert
        (
            gCelli,
            cellStencil(gCelli, move(neighbors), mesh_.cellCentres()[celli])
        );
        cellCells_[celli] = stencilMap_[gCelli];
    }

    if (Pstream::parRun())
    {
        Map<label> idxMap;
        List<labelList> sendMap(Pstream::nProcs());
        List<cellStencil> newCells;
        label idx = 0;
        forAll(cellCells_, i)
        {
            label celli = gIndex_.toGlobal(i);
            boolList added(Pstream::nProcs(), false);
            bool oneAdded = false;
            const cellStencil& lStencil(stencilMap_[celli]);

            forAll(lStencil, j)
            {
                label cellj = lStencil[j];
                label proci = gIndex_.whichProcID(cellj);
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

    //- Clean up map (remove non local indicies, and convert to local)
    forAllIter
    (
        Map<cellStencil>,
        stencilMap_,
        iter
    )
    {
        iter().updateLocalStencil(gIndex_);
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
            gIndex_,
            cellCells_,
            compactMap
        )
    );
}


// ************************************************************************* //
