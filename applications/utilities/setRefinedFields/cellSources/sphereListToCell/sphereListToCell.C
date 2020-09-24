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

#include "sphereListToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphereListToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, sphereListToCell, word);
    addToRunTimeSelectionTable(topoSetSource, sphereListToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::sphereListToCell::usage_
(
    sphereListToCell::typeName,
    "\n    Usage: sphereListToCell ((centreX centreY centreZ)) radii\n\n"
    "    Select all cells with cellCentre within bounding sphere\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sphereListToCell::combine
(
    topoSet& set,
    const vector& centre,
    const scalar& radius,
    const bool add
) const
{
    const pointField& ctrs = mesh_.cellCentres();

    const scalar radSquared = radius*radius;

    forAll(ctrs, celli)
    {
        scalar offset = magSqr(centre - ctrs[celli]);
        if (offset <= radSquared)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphereListToCell::sphereListToCell
(
    const polyMesh& mesh,
    const List<vector>& centres,
    const List<scalar>& radii
)
:
    topoSetSource(mesh),
    centres_(centres),
    radii_(radii)
{}


Foam::sphereListToCell::sphereListToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    centres_(dict.lookup("centres")),
    radii_
    (
        dict.found("radii")
      ? dict.lookup("radii")
      : List<scalar>(readScalar(dict.lookup("radius")), centres_.size())
    )
{}


Foam::sphereListToCell::sphereListToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    centres_(checkIs(is)),
    radii_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sphereListToCell::~sphereListToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sphereListToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    forAll(centres_, i)
    {
        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding cells with centre within sphere, with centre = "
                << centres_[i] << " and radius = " << radii_[i] << endl;

            combine(set, centres_[i], radii_[i], true);
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing cells with centre within sphere, with centre = "
                << centres_[i] << " and radius = " << radii_[i] << endl;

            combine(set, centres_[i], radii_[i], false);
        }
    }
}


// ************************************************************************* //
