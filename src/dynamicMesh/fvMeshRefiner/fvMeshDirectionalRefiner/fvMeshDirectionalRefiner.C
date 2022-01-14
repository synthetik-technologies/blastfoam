/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
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

#include "fvMeshDirectionalRefiner.H"
#include "newRefinementIterator.H"
#include "RefineBalanceMeshObject.H"
#include "parcelCloud.H"
#include "hexRef.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshDirectionalRefiner, 0);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshDirectionalRefiner, fvMesh);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshDirectionalRefiner, dictionary);
}


void Foam::fvMeshDirectionalRefiner::addCells
(
    const Map<label>& splitMap,
    List<refineCell>& refCells
)
{
    label newRefI = refCells.size();

    label oldSize = refCells.size();

    refCells.setSize(newRefI + splitMap.size());

    for (label refI = 0; refI < oldSize; refI++)
    {
        const refineCell& refCell = refCells[refI];

        Map<label>::const_iterator iter = splitMap.find(refCell.cellNo());

        if (iter == splitMap.end())
        {
            FatalErrorInFunction
                << "Problem : cannot find added cell for cell "
                << refCell.cellNo() << abort(FatalError);
        }

        refCells[newRefI++] = refineCell(iter(), refCell.direction());
    }
}


void Foam::fvMeshDirectionalRefiner::update
(
    const Map<label>& splitMap,
    vectorField& field
)
{
    field.setSize(field.size() + splitMap.size());

    forAllConstIter(Map<label>, splitMap, iter)
    {
        field[iter()] = field[iter.key()];
    }
}


void Foam::fvMeshDirectionalRefiner::addCells
(
    const Map<label>& splitMap,
    labelList& labels
)
{
    label newCelli = labels.size();

    labels.setSize(labels.size() + splitMap.size());

    forAllConstIter(Map<label>, splitMap, iter)
    {
        labels[newCelli++] = iter();
    }
}

Foam::vector Foam::fvMeshDirectionalRefiner::getDirection(const word& d)
{
    if (d == "x")
    {
        return vector(1, 0, 0);
    }
    else if (d == "y")
    {
        return vector(0, 1, 0);
    }
    else if (d == "z")
    {
        return vector(0, 0, 1);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown direction " << d << ". Valid direction are x, y, or z" << endl
            << abort(FatalError);
    }
    return vector::zero;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::fvMeshDirectionalRefiner::addCells
(
    const Map<label>& splitMap,
    labelListList& addedCells
)
{
    // Construct inverse addressing: from new to original cell.
    labelList origCell(mesh_.nCells(), -1);

    forAll(addedCells, celli)
    {
        const labelList& added = addedCells[celli];

        forAll(added, i)
        {
            label slave = added[i];

            if (origCell[slave] == -1)
            {
                origCell[slave] = celli;
            }
            else if (origCell[slave] != celli)
            {
                FatalErrorInFunction
                    << "Added cell " << slave << " has two different masters:"
                    << origCell[slave] << " , " << celli
                    << abort(FatalError);
            }
        }
    }


    forAllConstIter(Map<label>, splitMap, iter)
    {
        label masterI = iter.key();
        label newCelli = iter();

        while (origCell[masterI] != -1 && origCell[masterI] != masterI)
        {
            masterI = origCell[masterI];
        }

        if (masterI >= addedCells.size())
        {
            FatalErrorInFunction
                << "Map of added cells contains master cell " << masterI
                << " which is not a valid cell number" << endl
                << "This means that the mesh is not consistent with the"
                << " done refinement" << endl
                << "newCell:" << newCelli << abort(FatalError);
        }

        labelList& added = addedCells[masterI];

        if (added.empty())
        {
            added.setSize(2);
            added[0] = masterI;
            added[1] = newCelli;
        }
        else if (findIndex(added, newCelli) == -1)
        {
            label sz = added.size();
            added.setSize(sz + 1);
            added[sz] = newCelli;
        }
    }
}


Foam::labelListList Foam::fvMeshDirectionalRefiner::refineDirections
(
    List<vectorField>& cellDirections,
    List<labelList>& cellsToRefine
)
{
    labelListList addedCells(mesh_.nCells());

    // Iterator
    newRefinementIterator refiner(mesh_, cutter_(), cellWalker_());

    forAll(cellDirections, dirI)
    {
        if (debug)
        {
            Pout<< "multiDirRefinement : Refining " << cellsToRefine[dirI].size()
                << " cells in direction " << dirI << endl
                << endl;
        }

        const vectorField& dirField = cellDirections[dirI];

        // Combine cell to be cut with direction to cut in.
        // If dirField is only one element use this for all cells.

        List<refineCell> refCells(cellsToRefine[dirI].size());

        if (dirField.size() == 1)
        {
            // Uniform directions.
            if (debug)
            {
                Pout<< "fvMeshDirectionalRefiner : Uniform refinement:"
                    << dirField[0] << endl;
            }

            forAll(refCells, refI)
            {
                label celli = cellsToRefine[dirI][refI];

                refCells[refI] = refineCell(celli, dirField[0]);
            }
        }
        else
        {
            // Non uniform directions.
            forAll(refCells, refI)
            {
                const label celli = cellsToRefine[dirI][refI];

                refCells[refI] = refineCell(celli, dirField[refI]);
            }
        }

        // Do refine mesh (multiple iterations). Remember added cells.
        Map<label> splitMap = refiner.setRefinement(refCells);

        // Update overall map of added cells
        addCells(splitMap, addedCells);

        // Add added cells to list of cells to refine in next iteration
        addCells(splitMap, cellsToRefine[dirI]);

        // Update refinement direction for added cells.
        if (dirField.size() != 1)
        {
            forAll(cellDirections, i)
            {
                update(splitMap, cellDirections[i]);
            }
        }

        if (debug)
        {
            Pout<< "fvMeshDirectionalRefiner : Done refining direction " << dirI
                << " resulting in " << cellsToRefine[dirI].size() << " cells" << nl
                << endl;
        }
    }
    return addedCells;
}


Foam::labelListList
Foam::fvMeshDirectionalRefiner::refine(const labelList& origCellsToRefine)
{
    List<labelList> cellsToRefine;
    List<vectorField> cellDirections;
    if (refinementZones_.size())
    {
        labelHashSet selectedCells(origCellsToRefine);
        forAllConstIter
        (
            HashTable<List<vector>>,
            refinementZones_,
            iter
        )
        {
            label zoneID = mesh_.cellZones().findIndex(iter.key());
            if (zoneID == -1)
            {
                continue;
                FatalErrorInFunction
                    << "Could not find cellZone " << iter.key() << nl
                    << "cellZones:" << nl
                    << mesh_.cellZones().names() << endl
                    << abort(FatalError);
            }
            const cellZone& zone = mesh_.cellZones()[zoneID];
            DynamicList<label> cells(zone.size());
            forAll(zone, i)
            {
                if (selectedCells.found(zone[i]))
                {
                    cells.append(zone[i]);
                }
            }
            cells.shrink();
            const List<vector>& dirs(iter());

            forAll(dirs, diri)
            {
                cellsToRefine.append(cells);
                cellDirections.append(vectorField(1, dirs[diri]));
            }
        }
    }
    else
    {
        forAll(refinementDirections_, diri)
        {
            cellsToRefine.append(origCellsToRefine);
            cellDirections.append(vectorField(1, refinementDirections_[diri]));
        }
    }

    return refineDirections(cellDirections, cellsToRefine);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDirectionalRefiner::fvMeshDirectionalRefiner(fvMesh& mesh)
:
    fvMeshRefiner(mesh),

    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(mesh.nCells(), 0)
    ),

    cellWalker_(new geomCellLooper(mesh)),
    cutter_(new undoableMeshCutter(mesh, false)),

    isRefining_(false)
{
    Vector<label> geoD(mesh.geometricD());
    forAll(geoD, cmpti)
    {
        if (geoD[cmpti] > 0)
        {
            vector v(Zero);
            v[cmpti] = 1.0;
            refinementDirections_.append(v);
        }
    }
}


Foam::fvMeshDirectionalRefiner::fvMeshDirectionalRefiner
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
:
    fvMeshRefiner(mesh, dict, force, read),

    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh.time().timeName(),
            mesh,
            read ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        labelList(mesh.nCells(), 0)
    ),

    cellWalker_(),
    cutter_(new undoableMeshCutter(mesh, false)),

    isRefining_(false)
{
    // Read static part of dictionary
    readDict(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDirectionalRefiner::~fvMeshDirectionalRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshDirectionalRefiner::readDict(const dictionary& dict)
{
    fvMeshRefiner::readDict(dict);

    if (dict_.lookupOrDefault<Switch>("useGeoCut", false))
    {
        cellWalker_.reset(new geomCellLooper(mesh_));
    }
    else
    {
        cellWalker_.reset(new hexCellLooper(mesh_));
    }

    if (dict_.found("refinementZones"))
    {
        const PtrList<entry>& zones = dict_.lookup("refinementZones");
        forAll(zones, zonei)
        {
            const dictionary& refDict = zones[zonei].dict();
            refinementZones_.insert(zones[zonei].keyword(), {});

            List<vector>& dirs = refinementZones_[zones[zonei].keyword()];
            if (refDict.found("direction"))
            {
                dirs.append(getDirection(refDict.lookup<word>("direction")));
            }
            else if (refDict.found("directions"))
            {
                ITstream is(refDict.lookup("directions"));
                while (is.good())
                {
                    dirs.append(getDirection(word(is)));
                }
            }
        }
    }
    else if (dict_.found("direction"))
    {
        refinementDirections_.append(getDirection(dict_.lookup<word>("direction")));
    }
    else if (dict_.found("directions"))
    {
        ITstream is(dict_.lookup("directions"));
        while (is.good())
        {
            refinementDirections_.append(getDirection(word(is)));
        }
    }
    else
    {
        Vector<label> geoD(mesh_.geometricD());
        forAll(geoD, cmpti)
        {
            if (geoD[cmpti] > 0)
            {
                vector v(Zero);
                v[cmpti] = 1.0;
                refinementDirections_.append(v);
            }
        }
    }
}

bool Foam::fvMeshDirectionalRefiner::refine
(
    const scalarField& error,
    const labelList& maxCellLevel,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar unrefineLevel
)
{
    bool hasChanged = false;
    bool balanced = false;

    if (refineInterval_ == 0)
    {
        mesh_.topoChanging(hasChanged);

        return false;
    }

    if (mesh_.time().timeIndex() % refineInterval_ == 0)
    {
        HashTable<parcelCloud*> clouds
        (
            mesh_.lookupClass<parcelCloud>()
        );
        forAllIter(HashTable<parcelCloud*>, clouds, iter)
        {
            iter()->storeGlobalPositions();
        }

        labelList maxRefinement;
        if (!maxCellLevel.size())
        {
            maxRefinement.setSize
            (
                mesh_.nCells(),
                dict_.lookup<label>("maxRefinement")
            );
        }
        else
        {
            maxRefinement = maxCellLevel;
        }

        if (gMax(maxRefinement) == 0)
        {
            mesh_.topoChanging(hasChanged);

            return false;
        }
        else if (gMin(maxRefinement) < 0)
        {
            FatalErrorInFunction
                << "Illegal maximum refinement level " << gMin(maxRefinement) << nl
                << "The maxRefinement should be > 0." << nl
                << exit(FatalError);
        }

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(mesh_.nCells());

        // Determine candidates for refinement (looking at field only)
        selectRefineCandidates
        (
            lowerRefineLevel,
            upperRefineLevel,
            error,
            refineCell
        );

        // Extend with a buffer layer to prevent neighbouring points
        // being unrefined.
        for (label i = 0; i < nBufferLayers_; i++)
        {
            extendMarkedCells(refineCell, maxRefinement, i == 0);
        }

        if (mesh_.globalData().nTotalCells() < maxCells_)
        {
            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells_,
                    maxRefinement,
                    refineCell
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                labelListList oldToNew(refine(cellsToRefine));

                const labelList oldCellLevel(cellLevel_);
                cellLevel_.resize(mesh_.nCells());
                forAll(oldToNew, oldCelli)
                {
                    const labelList& added = oldToNew[oldCelli];

                    if (added.size())
                    {
                        forAll(added, i)
                        {
                            cellLevel_[added[i]] = oldCellLevel[oldCelli] + 1;
                        }
                    }
                    else
                    {
                        // Unrefined cell
                        cellLevel_[oldCelli] = oldCellLevel[oldCelli];
                    }
                }
                hasChanged = true;
            }
        }

        reduce(hasChanged, orOp<bool>());
        mesh_.topoChanging(hasChanged);
        balanced = balance();
        if (balanced)
        {
            //- Update objects stored on the mesh db
            BalanceMeshObject::updateObjects(mesh_);

            //- Update objects stores on the time db
            BalanceMeshObject::updateObjects
            (
                const_cast<Time&>(mesh_.time())
            );
            hasChanged = true;
        }
        else if (hasChanged)
        {
            //- Update objects stored on the mesh db
            RefineMeshObject::updateObjects(mesh_);

            //- Update objects stores on the time db
            RefineMeshObject::updateObjects
            (
                const_cast<Time&>(mesh_.time())
            );
        }

        if (hasChanged)
        {
            // Reset moving flag (if any). If not using inflation we'll not
            // move, if are using inflation any follow on movePoints will set
            // it.
            mesh_.moving(false);
        }
        nRefinementIterations_++;
    }

    return hasChanged || balanced;
}

void Foam::fvMeshDirectionalRefiner::updateMesh(const mapPolyMesh& map)
{
    fvMeshRefiner::updateMesh(map);
    if (!isRefining_)
    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexRef::updateMesh :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            hexRef::reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCelli)
            {
                label oldCelli = cellMap[newCelli];

                if (oldCelli == -1)
                {
                    newCellLevel[newCelli] = -1;
                }
                else
                {
                    newCellLevel[newCelli] = cellLevel_[oldCelli];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }
    }
}


void Foam::fvMeshDirectionalRefiner::distribute(const mapDistributePolyMesh& map)
{
    fvMeshRefiner::distribute(map);
    map.cellMap().distribute(cellLevel_);
}


bool Foam::fvMeshDirectionalRefiner::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool writeOk = cellLevel_.write();
    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            volScalarField::New
            (
                "cellLevel",
                mesh_,
                dimensionedScalar(dimless, 0)
            )
        );


        forAll(cellLevel_, celli)
        {
            scalarCellLevel[celli] = cellLevel_[celli];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// ************************************************************************* //
