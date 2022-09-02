/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "backupTopoSetSource.H"
#include "cellSet.H"
#include "cellZoneSet.H"
#include "faceSet.H"
#include "faceZoneSet.H"
#include "pointSet.H"
#include "pointZoneSet.H"
#include "coupledPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(backupTopoSetSource, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList
Foam::backupTopoSetSource::cellsToFaces(const labelList& cells) const
{
    PackedBoolList markedFaces(mesh_.nFaces());
    forAll(cells, ci)
    {
        const labelList& cFaces = mesh_.cells()[cells[ci]];
        forAll(cFaces, fi)
        {
            markedFaces.set(cFaces[fi]);
        }
    }
    syncTools::syncFaceList(mesh_, markedFaces, orEqOp<uint>());

    return markedFaces.used();
}


Foam::labelList
Foam::backupTopoSetSource::cellsToPoints(const labelList& cells) const
{
    PackedBoolList markedPoints(mesh_.nPoints());
    forAll(cells, ci)
    {
        const labelList& cPoints = mesh_.cellPoints()[cells[ci]];
        forAll(cPoints, pi)
        {
            markedPoints.set(cPoints[pi]);
        }
    }

    syncTools::syncPointList(mesh_, markedPoints, orEqOp<uint>(), 0);

    return markedPoints.used();
}



Foam::labelList
Foam::backupTopoSetSource::facesToCells(const labelList& faces) const
{
    if (!faces.size())
    {
        return labelList();
    }
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    labelHashSet selectedFaces(faces);
    labelHashSet selectedCells;
    forAll(owner, facei)
    {
        if (selectedFaces.found(owner[facei]))
        {
            if (facei >= mesh_.nInternalFaces())
            {
                selectedCells.insert(owner[facei]);
            }
            else if (selectedFaces.found(neighbour[facei]))
            {
                selectedCells.insert(owner[facei]);
            }
        }
    }
    return selectedCells.toc();
}


Foam::labelList
Foam::backupTopoSetSource::facesToPoints(const labelList& faces) const
{
    PackedBoolList markedPoints(mesh_.nPoints());
    forAll(faces, fi)
    {
        const labelList& fPoints = mesh_.faces()[faces[fi]];
        forAll(fPoints, pi)
        {
            markedPoints.set(fPoints[pi]);
        }
    }

    syncTools::syncPointList(mesh_, markedPoints, orEqOp<uint>(), false);
    return markedPoints.used();
}


Foam::labelList
Foam::backupTopoSetSource::pointsToCells(const labelList& points) const
{
    if (!points.size())
    {
        return labelList();
    }
    labelHashSet selectedPoints(points);
    labelHashSet selectedCells;
    forAll(mesh_.cells(), celli)
    {
        const cell& c = mesh_.cells()[celli];
        const labelList cp = c.labels(mesh_.faces());
        bool allFound = true;
        forAll(cp, pi)
        {
            if (!selectedPoints.found(cp[pi]))
            {
                allFound = false;
                break;
            }
        }
        if (allFound)
        {
            selectedCells.insert(celli);
        }
    }
    return selectedCells.toc();
}


Foam::labelList
Foam::backupTopoSetSource::pointsToFaces(const labelList& points) const
{
    if (!points.size())
    {
        return labelList();
    }
    labelHashSet selectedPoints(points);
    labelHashSet selectedFaces;
    forAll(mesh_.faces(), facei)
    {
        const face& f = mesh_.faces()[facei];
        bool allFound = true;
        forAll(f, pi)
        {
            if (!selectedPoints.found(f[pi]))
            {
                allFound = false;
                break;
            }
        }
        if (allFound)
        {
            selectedFaces.insert(facei);
        }
    }
    return selectedFaces.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backupTopoSetSource::backupTopoSetSource(const backupTopoSetSource& source)
:
    mesh_(source.mesh_),
    dict_(source.dict_),
    allowBackup_(source.allowBackup_),
    source_(source.source_->clone()),
    backup_
    (
        source.backup_.valid()
      ? source.backup_->clone()
      : autoPtr<topoSetSource>()
    )
{}


Foam::backupTopoSetSource::backupTopoSetSource
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    allowBackup_(dict.lookupOrDefault("allowBackup", false)),
    source_
    (
        topoSetSource::New
        (
            topoSetSourceType,
            mesh,
            dict
        )
    )
{
    if (dict.found("backup"))
    {
        backup_ =
            topoSetSource::New
            (
                topoSetSourceType,
                mesh,
                dict.subDict("backup")
            );
    }
}


Foam::backupTopoSetSource::backupTopoSetSource
(
    const polyMesh& mesh,
    const dictionary& dict,
    autoPtr<topoSetSource>& source
)
:
    mesh_(mesh),
    dict_(dict),
    allowBackup_(dict.lookupOrDefault("allowBackup", false)),
    source_(source),
    backup_()
{}


Foam::backupTopoSetSource::backupTopoSetSource
(
    const polyMesh& mesh,
    const dictionary& dict,
    autoPtr<topoSetSource>& source,
    autoPtr<topoSetSource>& backup
)
:
    mesh_(mesh),
    dict_(dict),
    allowBackup_(dict.lookupOrDefault("allowBackup", false)),
    source_(source),
    backup_(backup)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::backupTopoSetSource::~backupTopoSetSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::backupTopoSetSource::isCell() const
{
    return
        setType() == topoSetSource::CELLSETSOURCE
     || setType() == topoSetSource::CELLZONESOURCE;
}


bool Foam::backupTopoSetSource::isFace() const
{
    return
        setType() == topoSetSource::FACESETSOURCE
     || setType() == topoSetSource::FACEZONESOURCE;
}


bool Foam::backupTopoSetSource::isPoint() const
{
    return
        setType() == topoSetSource::POINTSETSOURCE
     || setType() == topoSetSource::POINTZONESOURCE;
}


void Foam::backupTopoSetSource::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    const label origSize = returnReduce(set.size(), sumOp<label>());
    source_->applyToSet(action, set);

    if
    (
        origSize == returnReduce(set.size(), sumOp<label>())
     && allowBackup_ && backup_.valid()
    )
    {
        backup_->applyToSet(action, set);
    }
}


void Foam::backupTopoSetSource::updateSets()
{
    autoPtr<topoSet> setPtr;
    switch (setType())
    {
        case topoSetSource::CELLSETSOURCE:
        {
            setPtr.set
            (
                new cellSet
                (
                    mesh_,
                    "cellSet",
                    mesh_.nCells()/10+1  // Reasonable size estimate.
                )
            );
            break;
        }
        case topoSetSource::CELLZONESOURCE:
        {
            setPtr.set
            (
                new cellSet
                (
                    mesh_,
                    "cellZoneSet",
                    mesh_.nCells()/10+1  // Reasonable size estimate.
                )
            );
            break;
        }
        case topoSetSource::FACESETSOURCE:
        {
            setPtr.set
            (
                new faceSet
                (
                    mesh_,
                    "faceSet",
                    mesh_.nFaces()/10+1  // Reasonable size estimate.
                )
            );
            break;
        }
        case topoSetSource::FACEZONESOURCE:
        {
            setPtr.set
            (
                new faceZoneSet
                (
                    mesh_,
                    "faceZoneSet",
                    mesh_.nFaces()/10+1  // Reasonable size estimate.
                )
            );
            break;
        }
        case topoSetSource::POINTSETSOURCE:
        {
            setPtr.set
            (
                new pointSet
                (
                    mesh_,
                    "pointSet",
                    mesh_.nPoints()/10+1  // Reasonable size estimate.
                )
            );
            break;
        }
        case topoSetSource::POINTZONESOURCE:
        {
            setPtr.set
            (
                new pointZoneSet
                (
                    mesh_,
                    "pointZoneSet",
                    mesh_.nPoints()/10+1  // Reasonable size estimate.
                )
            );
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown topoSetSource type " << setType() << endl
                << abort(FatalError);
        }
    }
    topoSet& set = setPtr();

    applyToSet(topoSetSource::NEW, set);

    if (isCell())
    {
        selectedCells_ = set.toc();
        selectedFaces_ = cellsToFaces(selectedCells_);
        selectedPoints_ = cellsToPoints(selectedCells_);
    }
    else if (isFace())
    {
//         if (setType() == topoSetSource::FACEZONESOURCE)
//         {
//             const faceZoneSet& fsz = dynamicCast<const faceZoneSet&>(set);
//             selectedFaces_ = fsz.addressing();
//             flipMap_ = fsz.flipMap();
//         }
//         else
        {
            selectedFaces_ = set.toc();
        }
        selectedCells_ = facesToCells(selectedFaces_);
        selectedPoints_ = facesToPoints(selectedFaces_);
    }
    else if (isPoint())
    {
        selectedPoints_ = set.toc();
        selectedCells_ = pointsToCells(selectedPoints_);
        selectedFaces_ = pointsToFaces(selectedPoints_);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type " << setType() << endl
            << abort(FatalError);
    }


    // Calculate the flip map
    if
    (
        returnReduce
        (
            flipMap_.size() != selectedFaces_.size(),
            orOp<bool>()
        )
    )
    {
        flipMap_ = boolList(selectedFaces_.size(), false);

        labelHashSet cells(selectedCells_);

        const labelList& owner = mesh_.faceOwner();
        const labelList& neighbour = mesh_.faceNeighbour();
        forAll(selectedFaces_, fi)
        {
            const label facei = selectedFaces_[fi];

            // Check id a face owner and neighbour are not the same for
            // internal faces
            if (facei < mesh_.nInternalFaces())
            {
                if
                (
                    !cells.found(owner[facei])
                 && cells.found(neighbour[facei])
                )
                {
                    flipMap_[fi] = true;
                }
            }
            // Assuming the faces have been synced, a face may be selected,
            // but not its owner cell
            else if (!cells.found(owner[facei]))
            {
                flipMap_[fi] = true;
            }
        }
    }
}


// ************************************************************************* //
