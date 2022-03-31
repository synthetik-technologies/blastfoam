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

#include "simpleCrackerFvMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleCohesiveZoneFvPatchVectorField.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleCrackerFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        simpleCrackerFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::simpleCrackerFvMesh::simpleCrackerFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleCrackerFvMesh::~simpleCrackerFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simpleCrackerFvMesh::update()
{
    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(*this);

    // Lookup displacement field
    const volVectorField* DPtr = NULL;
    if (solMod.incremental())
    {
        DPtr = &lookupObject<volVectorField>("DD");
    }
    else
    {
        DPtr = &lookupObject<volVectorField>("D");
    }
    const volVectorField& D = *DPtr;

    // Find cohesive patch and update crack

    label nFacesToBreak = 0;

    forAll(D.boundaryField(), patchI)
    {
        if
        (
            D.boundaryField()[patchI].type()
         == simpleCohesiveZoneFvPatchVectorField::typeName
        )
        {
            // Const-cast patch and call update crack
            simpleCohesiveZoneFvPatchVectorField& Dpatch =
                const_cast<simpleCohesiveZoneFvPatchVectorField&>
                (
                    refCast<const simpleCohesiveZoneFvPatchVectorField>
                    (
                        D.boundaryField()[patchI]
                    )
                );

            nFacesToBreak = Dpatch.updateCrack();

            break;
        }
    }

    Info<< nl << "Breaking " << nFacesToBreak << " faces" << nl << endl;

    if (time().outputTime())
    {
        polyMesh::write();
    }

    return bool(nFacesToBreak > 0);
}


// ************************************************************************* //
