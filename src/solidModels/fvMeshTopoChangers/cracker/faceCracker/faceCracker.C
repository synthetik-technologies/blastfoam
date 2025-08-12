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

Description
    Face cracker mesh modifier.  This modifier takes a set of
    internal face labels and converts them into boundary faces.

\*---------------------------------------------------------------------------*/

#include "faceCracker.H"
#include "polyTopoChange.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceCracker, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faceCracker::faceCracker
(
    const polyMesh& mesh,
    const word& faceZoneName,
    const word& crackPatchName
)
:
    mesh_(mesh),
    crackZone_(faceZoneName),
    crackPatch_(crackPatchName),
    coupledFacesToBreak_(),
    trigger_(false)
{
    if (mesh_.faceZones().findIndex(crackZone_) < 0)
    {
        FatalErrorInFunction
            << "Could not find crack zone " << crackZone_
            << ", valid zones are:" << nl
            << mesh_.faceZones() << endl
            << abort(FatalError);
    }
    if (mesh_.boundaryMesh().findIndex(crackPatch_) < 0)
    {
        FatalErrorInFunction
            << "Could not find crack patch " << crackPatch_
            << ", valid patches are:" << nl
            << mesh_.boundaryMesh().names() << endl
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "Face cracker object:" << nl
            << "    faceZone:   " << crackZone_ << nl
            << "    crackPatch: " << crackPatch_ << endl;
    }
}


// Construct from components
Foam::faceCracker::faceCracker
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    faceCracker
    (
        mesh,
        dict.lookup<word>("faceZone"),
        dict.lookup<word>("crackPatch")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceCracker::~faceCracker()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceCracker::setBreak
(
    polyMesh& mesh,
    const labelList& facesToBreak,
    const boolList& faceFlip,
    const labelList& coupledFacesToBreak
)
{
    if (trigger_)
    {
        FatalErrorInFunction
            << "Setting faces to break before previous break was "
            << "not executed.  Probably an error in topo change handling."
            << abort(FatalError);
    }

    // Check that all the faces in the face zone are internal
    if (debug)
    {
        // Check faces to break
        DynamicList<label> bouFacesInZone(facesToBreak.size());

        forAll(facesToBreak, faceI)
        {
            if (!mesh.isInternalFace(facesToBreak[faceI]))
            {
                bouFacesInZone.append(facesToBreak[faceI]);
            }
        }
    }

    // Put the faces into the face zone
    mesh.faceZones()[crackZone_].resetAddressing
    (
        facesToBreak,
        faceFlip
    );

    // Grab faces to break
    coupledFacesToBreak_ = coupledFacesToBreak;

    trigger_ = true;
}


bool Foam::faceCracker::changeTopology() const
{
    return trigger_;
}


void Foam::faceCracker::setRefinement(polyTopoChange& meshMod) const
{
    //detachFaceCracker(ref);

    detachInternalFaces(meshMod);
    detachCoupledFaces(meshMod);

    // Reset the trigger
    trigger_ = false;
}

// ************************************************************************* //
