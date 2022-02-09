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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(backupTopoSetSource, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList
Foam::backupTopoSetSource::cellsToFaces(const labelList& cells) const
{
    labelHashSet selectedFaces;
    forAll(cells, ci)
    {
        selectedFaces.insert(mesh_.cells()[cells[ci]]);
    }
    return selectedFaces.toc();
}


Foam::labelList
Foam::backupTopoSetSource::cellsToPoints(const labelList& cells) const
{
    labelHashSet selectedPoints;
    forAll(cells, ci)
    {
        selectedPoints.insert(mesh_.cellPoints()[cells[ci]]);
    }
    return selectedPoints.toc();
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
    labelHashSet selectedPoints;
    forAll(faces, fi)
    {
        selectedPoints.insert(mesh_.faces()[faces[fi]]);
    }
    return selectedPoints.toc();
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


void Foam::backupTopoSetSource::createSets
(
    labelList& cells,
    labelList& faces,
    boolList& flipMap,
    labelList& points
) const
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
        cells = set.toc();
        faces = cellsToFaces(cells);
        points = cellsToPoints(cells);
    }
    else if (isFace())
    {
        if (setType() == topoSetSource::FACEZONESOURCE)
        {
            const faceZoneSet& fsz = dynamicCast<const faceZoneSet&>(set);
            faces = fsz.addressing();
            flipMap = fsz.flipMap();
        }
        else
        {
            faces = set.toc();
        }
        cells = facesToCells(faces);
        points = facesToPoints(faces);
    }
    else if (isPoint())
    {
        points = set.toc();
        cells = pointsToCells(points);
        faces = pointsToFaces(points);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type " << setType() << endl
            << abort(FatalError);
    }
}


// ************************************************************************* //
