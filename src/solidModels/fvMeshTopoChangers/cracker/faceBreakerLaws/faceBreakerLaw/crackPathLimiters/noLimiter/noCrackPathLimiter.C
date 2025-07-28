/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "noCrackPathLimiter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace crackPathLimiters
{
    defineTypeNameAndDebug(none, 0);
    addToRunTimeSelectionTable(crackPathLimiter, none, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::crackPathLimiters::none::none
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    crackPathLimiter(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::crackPathLimiters::none::~none()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::crackPathLimiters::none::facesAllowedToBreak() const
{
    return surfaceScalarField::New
    (
        typeName + ":facesAllowedToBreak",
        mesh(),
        1.0
    );
}


void Foam::crackPathLimiters::none::clearOut()
{}


bool Foam::crackPathLimiters::none::write() const
{
    return true;
}

// ************************************************************************* //
