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

#include "hemisphericalMassToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hemisphericalMassToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, hemisphericalMassToCell, word);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hemisphericalMassToCell::hemisphericalMassToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    massToCell(dict),
    sphereToCell
    (
        mesh,
        centre(),
        R(dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hemisphericalMassToCell::~hemisphericalMassToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hemisphericalMassToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    labelHashSet oldSet(set);
    sphereToCell::applyToSet(action, set);

    checkMass
    (
        action == topoSetSource::REMOVE
      ? oldSet ^ set
      : (action == topoSetSource::ADD ? set ^ oldSet : set),
        mesh_
    );
}

// ************************************************************************* //
