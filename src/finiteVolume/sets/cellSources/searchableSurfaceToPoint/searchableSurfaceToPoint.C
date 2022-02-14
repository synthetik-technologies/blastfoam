/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "searchableSurfaceToPoint.H"
#include "polyMesh.H"
#include "searchableSurface.H"
#include "Time.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceToPoint, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToPoint,
        word
    );
}


Foam::topoSetSource::addToUsageTable Foam::searchableSurfaceToPoint::usage_
(
    searchableSurfaceToPoint::typeName,
    "\n    Usage: searchableSurfaceToPoint surface\n\n"
    "    Select all cells that are inside the surface "
    "\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceToPoint::combine
(
    topoSet& set,
    const bool add
) const
{
    List<volumeType> inOut(mesh_.nPoints(), volumeType::unknown);
    surfacePtr_->getVolumeType(mesh_.points(), inOut);
    forAll(inOut, pointi)
    {
        if (inOut[pointi] == volumeType::inside)
        {
            addOrDelete(set, pointi, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceToPoint::searchableSurfaceToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    surfacePtr_
    (
        searchableSurface::New
        (
            word(dict.lookup("surface")),
            IOobject
            (
                dict.lookupOrDefault("name", mesh.objectRegistry::db().name()),
                mesh.time().constant(),
                searchableSurface::geometryDir(mesh.time()),
                mesh.objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{
    if (!surfacePtr_->hasVolumeType())
    {
        FatalErrorInFunction
            << "Searchable surface type " << surfacePtr_->type()
            << " does not support volume type, but this is required" << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceToPoint::~searchableSurfaceToPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all points inside surface "
            << surfacePtr_().name() << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all points inside surface "
            << surfacePtr_().name() << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
