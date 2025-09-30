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

#include "triSurfaceRefPointToCell.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "Time.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceRefPointToCell, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        triSurfaceRefPointToCell,
        word
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurfaceRefPointToCell::combine
(
    topoSet& set,
    const bool add
) const
{
    List<pointIndexHit> info;
    tssPtr_->findNearest
    (
        mesh_.cellCentres(),
        scalarField(mesh_.nCells(), great),
        info
    );
    forAll(info, celli)
    {
        if (info[celli].hit())
        {
            const label facei = info[celli].index();
            const vector& n = surface_.faceNormals()[facei];
            const vector& fc = surface_.faceCentres()[facei];
            const vector& cc = mesh_.cellCentres()[celli];

            const vector diff(cc - fc);
            if (bb_.contains(cc) && (diff & n) > 0)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceRefPointToCell::triSurfaceRefPointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    surface_(fileName(dict.lookup("file"))),
    bb_(dict.lookupOrDefault("bounds", boundBox::greatBox)),
    tssPtr_(nullptr)
{
    const vector refPoint(dict.lookup("refPoint"));
    boolList includeFace(surface_.size());
    forAll(surface_, facei)
    {
        const vector& n = surface_.faceNormals()[facei];
        const vector& fc = surface_.faceCentres()[facei];
        const vector diff (refPoint - fc);
        const boundBox fBb(surface_.points(), surface_[facei], false);
        includeFace[facei] = bb_.overlaps(fBb) && (diff & n) > 0;
    }

    labelList faceMap, pointMap;
    surface_ = surface_.subsetMesh(includeFace, faceMap, pointMap);
    tssPtr_.set(new triSurfaceSearch(surface_));
    surface_.write("test.stl");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceRefPointToCell::~triSurfaceRefPointToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurfaceRefPointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells inside surface ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells inside surface ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
