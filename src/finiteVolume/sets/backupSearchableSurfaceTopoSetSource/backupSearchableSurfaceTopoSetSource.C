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

#include "backupSearchableSurfaceTopoSetSource.H"
#include "searchableSurfaceToCell.H"
#include "searchableSurfaceToPoint.H"
#include "cellSet.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(backupSearchableSurfaceTopoSetSource, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backupSearchableSurfaceTopoSetSource::
backupSearchableSurfaceTopoSetSource
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    backupTopoSetSource(topoSetSourceType, mesh, dict)
{
    if (isA<searchableSurfaceToCell>(source_()))
    {
        surfacePtr_.set
        (
            &dynamicCast<searchableSurfaceToCell>(source_()).surface()
        );
    }
    else if (isA<searchableSurfaceToPoint>(source_()))
    {
        surfacePtr_.set
        (
            &dynamicCast<searchableSurfaceToPoint>(source_()).surface()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Must use either "
            << searchableSurfaceToCell::typeName << " or "
            << searchableSurfaceToPoint::typeName << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::backupSearchableSurfaceTopoSetSource::~backupSearchableSurfaceTopoSetSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
