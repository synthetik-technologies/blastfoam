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

#include "backupSearchableSurface.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(backupSearchableSurface, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backupSearchableSurface::backupSearchableSurface
(
    const backupSearchableSurface& surface
)
:
    mesh_(surface.mesh_),
    dict_(surface.dict_),
    allowBackup_(false),
    surfacePtr_(surface.surfacePtr_->clone()),
    backupPtr_
    (
        surface.backupPtr_.valid()
      ? surface.backupPtr_->clone()
      : autoPtr<searchableSurface>()
    )
{}


Foam::backupSearchableSurface::backupSearchableSurface
(
    const word& surfaceType,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    allowBackup_(false),
    surfacePtr_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                surfaceType,
                mesh.time().constant(),
                mesh
            ),
            dict
        )
    ),
    backupPtr_()
{
    if (dict.isDict("backup"))
    {
        const dictionary& backupDict(dict.subDict("backup"));
        if (backupDict.found("surface"))
        {
            backupPtr_ =
                searchableSurface::New
                (
                    backupDict.lookup<word>("surface"),
                    IOobject
                    (
                        surfaceType + "_backup",
                        mesh.time().constant(),
                        mesh
                    ),
                    backupDict
                );
        }
    }
}


Foam::backupSearchableSurface::backupSearchableSurface
(
    searchableSurface* surface,
    const dictionary& dict
)
:
    mesh_(dynamicCast<const polyMesh>(surface->db())),
    dict_(dict),
    allowBackup_(false),
    surfacePtr_(surface),
    backupPtr_()
{
    if (dict.isDict("backup"))
    {
        const dictionary& backupDict(dict.subDict("backup"));
        if (backupDict.found("surface"))
        {
            backupPtr_ =
                searchableSurface::New
                (
                    backupDict.lookup<word>("surface"),
                    IOobject
                    (
                        surface->name() + "_backup",
                        mesh_.time().constant(),
                        mesh_
                    ),
                    backupDict
                );
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::backupSearchableSurface::~backupSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::backupSearchableSurface::selectPoints
(
    const pointField& points,
    const bool backupAllowed
) const
{
    DynamicList<label> selected(points.size());
    List<volumeType> inOut(points.size(), volumeType::unknown);
    surfacePtr_->getVolumeType(points, inOut);
    forAll(inOut, pointi)
    {
        if (inOut[pointi] == volumeType::inside)
        {
            selected.append(pointi);
        }
    }
    if
    (
        !returnReduce(selected.size(), sumOp<label>())
     && allowBackup_
     && backupAllowed
     && backupPtr_.valid()
    )
    {
        backupPtr_->getVolumeType(points, inOut);
        forAll(inOut, pointi)
        {
            if (inOut[pointi] == volumeType::inside)
            {
                selected.append(pointi);
            }
        }
    }
    return move(selected);
}


// ************************************************************************* //
