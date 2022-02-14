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

#include "searchableSurfaceToCell.H"
#include "polyMesh.H"
#include "searchableSurface.H"
#include "Time.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceToCell, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToCell,
        word
    );
}


Foam::topoSetSource::addToUsageTable Foam::searchableSurfaceToCell::usage_
(
    searchableSurfaceToCell::typeName,
    "\n    Usage: searchableSurfaceToCell surface\n\n"
    "    Select all cells that are inside the surface "
    "\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceToCell::combine
(
    topoSet& set,
    const bool add
) const
{
    List<volumeType> inOut(mesh_.nCells(), volumeType::unknown);
    surfacePtr_->getVolumeType(mesh_.cellCentres(), inOut);
    forAll(inOut, celli)
    {
        if (inOut[celli] == volumeType::inside)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceToCell::searchableSurfaceToCell
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

Foam::searchableSurfaceToCell::~searchableSurfaceToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells inside surface "
            << surfacePtr_().name() << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells inside surface "
            << surfacePtr_().name() << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
