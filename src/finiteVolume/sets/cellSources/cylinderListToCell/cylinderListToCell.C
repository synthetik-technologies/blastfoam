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

#include "cylinderListToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderListToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, cylinderListToCell, word);
    addToRunTimeSelectionTable(topoSetSource, cylinderListToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cylinderListToCell::usage_
(
    cylinderListToCell::typeName,
    "\n    Usage: cylinderListToCell ((p1X p1Y p1Z)) ((p2X p2Y p2Z)) radii\n\n"
    "    Select all cells with cell centres within bounding cylinders\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cylinderListToCell::combine
(
    topoSet& set,
    const vector& p1,
    const vector& p2,
    const scalar& radius,
    const bool add
) const
{
    const pointField& ctrs = mesh_.cellCentres();

    const vector axis = p2 - p1;
    const scalar rad2 = sqr(radius);
    const scalar magAxis2 = magSqr(axis);

    forAll(ctrs, celli)
    {
        vector d = ctrs[celli] - p1;
        scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxis2))
        {
            scalar d2 = (d & d) - sqr(magD)/magAxis2;
            if (d2 < rad2)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderListToCell::cylinderListToCell
(
    const polyMesh& mesh,
    const List<vector>& p1s,
    const List<vector>& p2s,
    const List<scalar> radii
)
:
    topoSetSource(mesh),
    p1s_(p1s),
    p2s_(p2s),
    radii_(radii)
{}


Foam::cylinderListToCell::cylinderListToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    p1s_(dict.lookup("p1s")),
    p2s_(dict.lookup("p2s")),
    radii_
    (
        dict.found("radii")
      ? dict.lookup("radii")
      : List<scalar>(readScalar(dict.lookup("radius")), p1s_.size())
    )
{}


Foam::cylinderListToCell::cylinderListToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    p1s_(checkIs(is)),
    p2s_(checkIs(is)),
    radii_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cylinderListToCell::~cylinderListToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylinderListToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    forAll(p1s_, i)
    {
        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding cells with centre within cylinder, with p1 = "
                << p1s_[i] << ", p2 = " << p2s_[i] << " and radius = "
                << radii_[i] << endl;

            combine(set, p1s_[i], p2s_[i], radii_[i], true);
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing cells with centre within cylinder, with p1 = "
                << p1s_[i] << ", p2 = " << p2s_[i] << " and radius = "
                << radii_[i] << endl;

            combine(set, p1s_[i], p2s_[i], radii_[i], false);
        }
    }
}


// ************************************************************************* //
