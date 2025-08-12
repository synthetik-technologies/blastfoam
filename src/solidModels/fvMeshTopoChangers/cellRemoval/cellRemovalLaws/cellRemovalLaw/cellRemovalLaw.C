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

Class
    cellRemovalLaw

\*---------------------------------------------------------------------------*/

#include "cellRemovalLaw.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellRemovalLaw, 0);
    defineRunTimeSelectionTable(cellRemovalLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellRemovalLaw::cellRemovalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    exposedPatch_(dict.lookup("exposedPatch"))
{
    if (exposedFacesPatchID() < 0)
    {
        FatalIOErrorInFunction(dict)
            << exposedPatch_ << " is not a valid patch. Valid patches are " << nl
            << mesh.boundaryMesh().names() << endl
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::cellRemovalLaw::~cellRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::cellRemovalLaw::exposedFacesPatchID() const
{
    return mesh_.boundaryMesh()[exposedPatch_].index();
}


// ************************************************************************* //
