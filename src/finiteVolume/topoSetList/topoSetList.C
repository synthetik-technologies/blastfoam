/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "topoSetList.H"
#include "topoSetSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSetList, 0);

    template<>
    const char* NamedEnum<topoSetList::SelectionType, 5>::names[] =
    {
        "all",
        "internal",
        "interface",
        "boundary",
        "interfaceAndBoundary"
    };

    const NamedEnum<topoSetList::SelectionType, 5>
        topoSetList::selectionNames;
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::topoSetList::topoSetList(const fvMesh& mesh)
:
    TopoSetList
    (
        mesh,
        IOobject
        (
            typeName,
            mesh.facesInstance(),
            mesh
        )
    )
{
    wordList cellZoneNames = mesh.cellZones().names();
    forAll(mesh.cellZones(), czi)
    {
        cellTopoSets_.insert
        (
            cellZoneNames[czi],
            new cellSet
            (
                mesh_,
                cellZoneNames[czi],
                labelHashSet(mesh.cellZones()[czi])
            )
        );
        cellZones_.insert(cellZoneNames[czi]);
    }

    wordList faceZoneNames = mesh.faceZones().names();
    forAll(mesh.faceZones(), fzi)
    {
        faceTopoSets_.insert
        (
            faceZoneNames[fzi],
            new faceSet
            (
                mesh_,
                faceZoneNames[fzi],
                labelHashSet(mesh.faceZones()[fzi])
            )
        );
        faceZones_.insert(faceZoneNames[fzi]);
    }

    wordList pointZoneNames = mesh.pointZones().names();
    forAll(mesh.pointZones(), pzi)
    {
        pointTopoSets_.insert
        (
            pointZoneNames[pzi],
            new pointSet
            (
                mesh_,
                pointZoneNames[pzi],
                labelHashSet(mesh.pointZones()[pzi])
            )
        );
        pointZones_.insert(pointZoneNames[pzi]);
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::topoSetList::~topoSetList()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

void Foam::topoSetList::updateCells
(
    const dictionary& dict,
    const labelList& selectedCells,
    const bool isZone
)
{
    word setType = "cell" + word(isZone ? "Zones" : "Sets");
    if (!dict.found(setType))
    {
        return;
    }
    PtrList<dictionary> setDicts(dict.lookup(setType));
    forAll(setDicts, seti)
    {
        modifyTopoSet
        (
            "cell",
            isZone,
            cellTopoSets_,
            cellZones_,
            cellSets_,
            selectedCells,
            setDicts[seti]
        );
    }
}


void Foam::topoSetList::updateFaces
(
    const dictionary& dict,
    const labelList& selectedFaces,
    const bool isZone
)
{
    word setType = "face" + word(isZone ? "Zones" : "Sets");
    if (!dict.found(setType))
    {
        return;
    }
    PtrList<dictionary> setDicts(dict.lookup(setType));
    forAll(setDicts, seti)
    {
        modifyTopoSet
        (
            "face",
            isZone,
            faceTopoSets_,
            faceZones_,
            faceSets_,
            extractSelectedFaces
            (
                setDicts[seti],
                selectedFaces,
                false
            ),
            setDicts[seti]
        );
    }
}


void Foam::topoSetList::updatePoints
(
    const dictionary& dict,
    const labelList& selectedPoints,
    const bool isZone
)
{
    word setType = "point" + word(isZone ? "Zones" : "Sets");
    if (!dict.found(setType))
    {
        return;
    }
    PtrList<dictionary> setDicts(dict.lookup(setType));
    forAll(setDicts, seti)
    {
        modifyTopoSet
        (
            "point",
            isZone,
            pointTopoSets_,
            pointZones_,
            pointSets_,
            extractSelectedPoints
            (
                setDicts[seti],
                selectedPoints,
                false
            ),
            setDicts[seti]
        );
    }
}


void Foam::topoSetList::update
(
    const dictionary& dict,
    const labelList& selectedCells,
    const labelList& selectedFaces,
    const labelList& selectedPoints
)
{
    // Sets
    updateCells(dict, selectedCells, false);
    updateFaces(dict, selectedFaces, false);
    updatePoints(dict, selectedPoints, false);

    // Zones
    updateCells(dict, selectedCells, true);
    updateFaces(dict, selectedFaces, true);
    updatePoints(dict, selectedPoints, true);
}


void Foam::topoSetList::modifyTopoSet
(
    const word& type,
    const bool isZone,
    HashPtrTable<topoSet>& topoSets,
    wordHashSet& zones,
    wordHashSet& sets,
    const labelList& selected,
    const dictionary& dict
)
{
    const word setName(dict.lookup("name"));
    const bool isSet(!isZone);
    const word setType = type + (isZone ? "ZoneSet" : "Set");

    topoSetSource::setAction action = topoSetSource::toAction
    (
        dict.lookup<word>("action")
    );

    if (action == topoSetSource::NEW)
    {
        topoSets.set
        (
            setName,
            topoSet::New(setType, mesh_, setName, selected.size()).ptr()
        );
    }
    else if (action == topoSetSource::CLEAR)
    {
        topoSets.set
        (
            setName,
            topoSet::New(setType, mesh_, setName, 0).ptr()
        );
    }
    else
    {
        if (!topoSets.found(setName))
        {
            if (dict.found("source"))
            {
                const word sourceName(dict.lookup("source"));
                if (!topoSets.found(sourceName))
                {
                    topoSets.set
                    (
                        sourceName,
                        topoSet::New
                        (
                            setType,
                            mesh_,
                            sourceName,
                            IOobject::MUST_READ
                        ).ptr()
                    );
                }
                topoSets.set
                (
                    setName,
                    topoSet::New
                    (
                        setType,
                        mesh_,
                        setName,
                        *topoSets[sourceName]
                    ).ptr()
                );
            }
            else
            {
                Info<< "source was not specified so reading " << setName
                    << endl;
                topoSets.set
                (
                    setName,
                    topoSet::New
                    (
                        setType,
                        mesh_,
                        setName,
                        IOobject::MUST_READ
                    ).ptr()
                );
                sets.insert(setName);
            }
        }
    }
    topoSet& currentSet = *topoSets[setName];

    // Handle special actions (clear, invert) locally, rest through
    // sources.
    switch (action)
    {
        case topoSetSource::NEW:
        case topoSetSource::ADD:
        {
            forAll(selected, i)
            {
                currentSet.insert(selected[i]);
            }
            if (isZone)
            {
                zones.insert(setName);
            }
            else if (isSet)
            {
                sets.insert(setName);
            }
        }
        break;

        case topoSetSource::DELETE:
        {
            HashPtrTable<topoSet>::iterator iter = topoSets.find(setName);
            topoSets.erase(iter);
        }
        break;

        case topoSetSource::SUBSET:
        {
            // Backup current set.
            autoPtr<topoSet> oldSet
            (
                topoSet::New
                (
                    setType,
                    mesh_,
                    currentSet.name() + "_old2",
                    selected.size()
                )
            );
            forAll(selected, i)
            {
                oldSet().insert(selected[i]);
            }

            const word sourceName(dict.lookup("source"));
            const topoSet* sourcePtr = nullptr;
            if (!topoSets.found(sourceName))
            {
                sourcePtr = topoSet::New
                (
                    setType,
                    mesh_,
                    sourceName,
                    IOobject::MUST_READ
                ).ptr();
            }
            else
            {
                sourcePtr = topoSets[sourceName];
            }
            currentSet = *sourcePtr;

            if (!topoSets.found(sourceName))
            {
                deleteDemandDrivenData(sourcePtr);
            }

            // Combine new value of currentSet with old one.
            currentSet.subset(oldSet);

            // Synchronise for coupled patches.
            currentSet.sync(mesh_);

            if (isZone)
            {
                zones.insert(setName);
            }
            else if (isSet)
            {
                sets.insert(setName);
            }
        }
        break;

        case topoSetSource::CLEAR:
        break;

        case topoSetSource::INVERT:
            topoSets[setName]->invert(topoSets[setName]->maxSize(mesh_));
            if (isZone)
            {
                zones.insert(setName);
            }
            else if (isSet)
            {
                sets.insert(setName);
            }
        break;

        case topoSetSource::REMOVE:
        {
            forAll(selected, i)
            {
                currentSet.erase(selected[i]);
            }

            if (isZone)
            {
                zones.insert(setName);
            }
            else if (isSet)
            {
                sets.insert(setName);
            }
        }
        break;


        default:
            WarningInFunction
                << "Unhandled action " << action << endl;
        break;
    }
}


Foam::labelList Foam::topoSetList::extractInterfaceCells
(
    const labelList& cells
) const
{
    labelHashSet cellSet(cells);
    labelList newCells(cells);

    label I = 0;
    forAll(cells, ci)
    {
        const label celli = cells[ci];
        const labelList& cc = mesh_.cellCells()[celli];
        bool found = cellSet.found(celli);
        bool allSame = true;
        forAll(cc, cj)
        {
            if (cellSet.found(cc[cj]) != found)
            {
                allSame = false;
                break;
            }
        }
        if (!allSame)
        {
            newCells[I++] = celli;
        }
    }
    newCells.resize(I);
    return newCells;
}


//- Remove faces without an owner and a neighbour
Foam::labelList Foam::topoSetList::extractSelectedFaces
(
    const dictionary& dict,
    const labelList& faces,
    const bool defaultAll
) const
{
    const SelectionType sType =
        defaultAll
      ? selectionNames[dict.lookupOrDefault<word>("selectionMode", "all")]
      : selectionNames[dict.lookup<word>("selectionMode")];

    // Return if the selection is empty of all faces are selected
    if (!returnReduce(faces.size(), sumOp<label>()) || sType == ALL)
    {
        return faces;
    }

    // Create an easily searchable list of selected faces for searching
    const labelHashSet selectedFaces(faces);

    // Add internal faces (coupled faces are included)
    labelHashSet internalFaces(faces.size());
    if (sType == INTERNAL)
    {
        forAll(faces, fi)
        {
            const label facei = faces[fi];
            if (facei >= mesh_.nInternalFaces())
            {
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);
                if (mesh_.boundaryMesh()[patchi].coupled())
                {
                    internalFaces.insert(faces[fi]);
                }
            }
            else
            {
                internalFaces.insert(faces[fi]);
            }
        }
        return internalFaces.toc();
    }

    // Add boundary faces on the selected faces
    labelHashSet boundaryFaces(faces.size());
    if (sType == BOUNDARY || sType == INTERFACE_AND_BOUNDARY)
    {
        wordReList patchNames(dict.lookup("patches"));

        forAll(patchNames, i)
        {
            const label patchi =
                mesh_.boundaryMesh().findIndex(patchNames[i]);
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];

            forAll(pp, fi)
            {
                const label facei = pp.start() + fi;
                if (selectedFaces.found(facei))
                {
                    boundaryFaces.insert(facei);
                }
            }
        }
        if (sType == BOUNDARY)
        {
            return boundaryFaces.toc();
        }
    }

    // Select interface faces, i.e. selected faces whos owner/neighbour
    // cells have some faces that are not selected
    labelHashSet interfaceFaces(faces.size());
    if (sType == INTERFACE || sType == INTERFACE_AND_BOUNDARY)
    {
        const labelList& owner = mesh_.faceOwner();
        const labelList& neighbour = mesh_.faceNeighbour();
        forAll(faces, fi)
        {
            const label facei = faces[fi];
            const cell& own = mesh_.cells()[owner[facei]];

            bool added = false;
            forAll(own, fj)
            {
                if (!selectedFaces.found(own[fj]))
                {
                    interfaceFaces.insert(facei);
                    break;
                }
            }
            if (!added && facei < mesh_.nInternalFaces())
            {
                const cell& nei = mesh_.cells()[neighbour[facei]];
                forAll(nei, fj)
                {
                    if (!selectedFaces.found(nei[fj]))
                    {
                        interfaceFaces.insert(facei);
                        break;
                    }
                }
            }
        }

        if (sType == INTERFACE)
        {
            return interfaceFaces.toc();
        }
    }

    return (interfaceFaces | boundaryFaces).toc();
}


Foam::labelList Foam::topoSetList::extractSelectedPoints
(
    const dictionary& dict,
    const labelList& points,
    const bool defaultAll
) const
{
    const SelectionType sType =
        defaultAll
      ? selectionNames[dict.lookupOrDefault<word>("selectionMode", "all")]
      : selectionNames[dict.lookup("selectionMode")];

    // Return if the selection is empty of all points are selected
    if (!returnReduce(points.size(), sumOp<label>()) || sType == ALL)
    {
        return points;
    }

    // Create a searchable list of selected points
    const labelHashSet selectedPoints(points);

    // Select all boundary points
    labelHashSet boundaryPoints;
    if (sType == BOUNDARY || sType == INTERFACE_AND_BOUNDARY)
    {
        // Only select points on the given patches
        labelHashSet patchIDs
        (
            mesh_.boundaryMesh().patchSet
            (
                dict.lookup<wordReList>("patches")
            )
        );

        // Loop through all the selected patches and add the selected points
        forAll(patchIDs, i)
        {
            const label patchi = patchIDs[i];
            const labelList& ppoints
            (
                mesh_.boundaryMesh()[patchi].meshPoints()
            );
            forAll(ppoints, pi)
            {
                if (selectedPoints.found(ppoints[pi]))
                {
                    boundaryPoints.insert(ppoints[pi]);
                }
            }
        }

        if (sType == BOUNDARY)
        {
            return boundaryPoints.toc();
        }
    }

    // The internal points are the entire selection minus the boundary points
    if (sType == INTERNAL)
    {
        labelHashSet internalPoints(points);
        internalPoints -= boundaryPoints;
        return internalPoints.toc();
    }


    // Select the interface points, i.e. selected points for neighbour
    // points that are not selected
    const labelListList& pointPoints = mesh_.pointPoints();
    labelHashSet interfacePoints;
    forAll(points, pi)
    {
        const label pointi = points[pi];
        const labelList& pp = pointPoints[pointi];

        bool add = false;
        forAll(pp, pj)
        {
            if (!selectedPoints.found(pp[pj]))
            {
                add = true;
                break;
            }
        }
        if (add)
        {
            interfacePoints.insert(pointi);
        }
    }

    if (sType == INTERFACE)
    {
        return interfacePoints.toc();
    }
    else if (sType == INTERFACE_AND_BOUNDARY)
    {
        return (boundaryPoints | interfacePoints).toc();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown selection mode " << dict.lookup<word>("selectionMode")
            << ", valid selection modes are: " << nl
            << selectionNames.toc() << endl
            << abort(FatalError);
    }
    return points;
}


void Foam::topoSetList::clear()
{
    cellTopoSets_.clear();
    faceTopoSets_.clear();
    pointTopoSets_.clear();

    cellSets_.clear();
    cellZones_.clear();

    faceSets_.clear();
    faceZones_.clear();

    pointSets_.clear();
    faceZones_.clear();
}



void Foam::topoSetList::updateMesh(const mapPolyMesh& morphMap)
{
    forAllConstIter(HashPtrTable<topoSet>, cellTopoSets_, iter)
    {
        iter()->updateMesh(morphMap);
    }
    forAllConstIter(HashPtrTable<topoSet>, faceTopoSets_, iter)
    {
        iter()->updateMesh(morphMap);
    }
    forAllConstIter(HashPtrTable<topoSet>, pointTopoSets_, iter)
    {
        iter()->updateMesh(morphMap);
    }
}


void Foam::topoSetList::transferZones(const bool remove)
{
    List<cellZone*> meshCellZones(cellZones_.size());
    List<faceZone*> meshFaceZones(faceZones_.size());
    List<pointZone*> meshPointZones(pointZones_.size());

    wordList cellZones(cellZones_.toc());
    forAll(cellZones, zonei)
    {
        const word& zoneName = cellZones[zonei];
        HashPtrTable<topoSet>::iterator czIter = cellTopoSets_.find(zoneName);
        if (czIter != cellTopoSets_.end())
        {
            meshCellZones[zonei] =
                new cellZone
                (
                    czIter()->name(),
                    czIter()->toc(),
                    zonei,
                    mesh_.cellZones()
                );
            if (remove)
            {
                cellTopoSets_.erase(czIter);
            }

        }
    }

    wordList faceZones(faceZones_.toc());
    forAll(faceZones, zonei)
    {
        const word& zoneName = faceZones[zonei];
        HashPtrTable<topoSet>::iterator fzIter = faceTopoSets_.find(zoneName);
        if (fzIter != faceTopoSets_.end())
        {
            meshFaceZones[zonei] =
                new faceZone
                (
                    fzIter()->name(),
                    fzIter()->toc(),
                    boolList(fzIter()->size(), false),
                    zonei,
                    mesh_.faceZones()
                );
            if (remove)
            {
                faceTopoSets_.erase(fzIter);
            }

        }
    }

    wordList pointZones(pointZones_.toc());
    forAll(pointZones, zonei)
    {
        const word& zoneName = pointZones[zonei];
        HashPtrTable<topoSet>::iterator pzIter = pointTopoSets_.find(zoneName);
        if (pzIter != pointTopoSets_.end())
        {
            meshPointZones[zonei] =
                new pointZone
                (
                    pzIter()->name(),
                    pzIter()->toc(),
                    zonei,
                    mesh_.pointZones()
                );
            if (remove)
            {
                pointTopoSets_.erase(pzIter);
            }

        }
    }
    if (remove)
    {
        cellZones_.clear();
        faceZones_.clear();
        pointZones_.clear();
    }


    if (meshCellZones.size() || meshFaceZones.size() || meshPointZones.size())
    {
        if (meshCellZones.size() && debug)
        {
            Info << "Adding cellZones " << cellZones_ << endl;
        }
        if (meshFaceZones.size() && debug)
        {
            Info << "Adding faceZones " << faceZones_ << endl;
        }
        if (meshPointZones.size() && debug)
        {
            Info << "Adding pointZones " << pointZones_ << endl;
        }
        fvMesh& mesh = const_cast<fvMesh&>(mesh_);
        mesh.pointZones().clear();
        mesh.faceZones().clear();
        mesh.cellZones().clear();

        mesh.addZones(meshPointZones, meshFaceZones, meshCellZones);
    }
}


bool Foam::topoSetList::writeSets() const
{
    List<cellZone*> meshCellZones;
    List<faceZone*> meshFaceZones;
    List<pointZone*> meshPointZones;
    label cellZonei = 0;
    forAllConstIter
    (
        HashPtrTable<topoSet>,
        cellTopoSets_,
        iter
    )
    {
        if (cellZones_.found(iter()->name()))
        {
            meshCellZones.append
            (
                new cellZone
                (
                    iter()->name(),
                    iter()->toc(),
                    cellZonei++,
                    mesh_.cellZones()
                )
            );
        }
        else if (cellSets_.found(iter()->name()))
        {
            DebugInfo<< "Writing cell set " << iter()->name() << endl;
            iter()->instance() = mesh_.polyMesh::instance();
            iter()->write();
        }
    }

    label faceZonei = 0;
    forAllConstIter
    (
        HashPtrTable<topoSet>,
        faceTopoSets_,
        iter
    )
    {
        if (faceZones_.found(iter()->name()))
        {
            meshFaceZones.append
            (
                new faceZone
                (
                    iter()->name(),
                    iter()->toc(),
                    boolList(iter()->size(), false),
                    faceZonei++,
                    mesh_.faceZones()
                )
            );
        }
        else if (faceSets_.found(iter()->name()))
        {
            DebugInfo<< "Writing face set " << iter()->name() << endl;
            iter()->instance() = mesh_.polyMesh::instance();
            iter()->write();
        }
    }

    label pointZonei = 0;
    forAllConstIter
    (
        HashPtrTable<topoSet>,
        pointTopoSets_,
        iter
    )
    {
        if (pointZones_.found(iter()->name()))
        {
            meshPointZones.append
            (
                new pointZone
                (
                    iter()->name(),
                    iter()->toc(),
                    pointZonei++,
                    mesh_.pointZones()
                )
            );
        }
        else if (pointSets_.found(iter()->name()))
        {
            DebugInfo<< "Writing point set " << iter()->name() << endl;
            iter()->instance() = mesh_.polyMesh::instance();
            iter()->write();
        }
    }
    if (cellZonei || faceZonei || pointZonei)
    {
        if (cellZonei && debug)
        {
            Info << "Adding cellZones " << cellZones_ << endl;
        }
        if (faceZonei && debug)
        {
            Info << "Adding faceZones " << faceZones_ << endl;
        }
        if (pointZonei && debug)
        {
            Info << "Adding pointZones " << pointZones_ << endl;
        }
        fvMesh& mesh = const_cast<fvMesh&>(mesh_);
        mesh.pointZones().clear();
        mesh.faceZones().clear();
        mesh.cellZones().clear();

        mesh.addZones(meshPointZones, meshFaceZones, meshCellZones);
        return true;
    }
    return false;
}


// ************************************************************************* //
