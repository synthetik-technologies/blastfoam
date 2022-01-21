/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Virtual base class for mesh modifiers.

\*---------------------------------------------------------------------------*/

#include "polyMeshModifier.H"
#include "dictionary.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshModifier, 0);

    defineRunTimeSelectionTable(polyMeshModifier, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyMeshModifier::polyMeshModifier
(
    const word& name,
    const label index,
    const polyTopoChanger& mme,
    const bool act
)
:
    name_(name),
    index_(index),
    topoChanger_(mme),
    active_(act)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMeshModifier::~polyMeshModifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyTopoChanger& Foam::polyMeshModifier::topoChanger() const
{
    return topoChanger_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyMeshModifier& pmm)
{
    pmm.write(os);
    os.check("Ostream& operator<<(Ostream& f, const polyMeshModifier& pmm)");
    return os;
}


// ************************************************************************* //
