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
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceCracker, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        faceCracker,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceCracker::checkDefinition()
{
    if
    (
        !crackZoneID_.active()
     || !crackPatchID_.active()
    )
    {
        FatalErrorIn
        (
            "void Foam::faceCracker::checkDefinition()"
        )   << "Not all zones and patches needed in the definition "
            << "have been found.  Please check your mesh definition.\n"
            << "\tcrackZoneID_.active() is " << crackZoneID_.active() << nl
            << "\tcrackPatchID_.active() is " << crackPatchID_.active() //<< nl
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "Face cracker object " << name() << " :" << nl
            << "    faceZoneID:   " << crackZoneID_ << nl
            << "    crackPatchID: " << crackPatchID_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faceCracker::faceCracker
(
    const word& name,
    const label index,
    const polyTopoChanger& mme,
    const word& faceZoneName,
    const word& crackPatchName
)
:
    polyMeshModifier(name, index, mme, true),
    crackZoneID_(faceZoneName, mme.mesh().faceZones()),
    crackPatchID_(crackPatchName, mme.mesh().boundaryMesh()),
    coupledFacesToBreak_(),
    trigger_(false)
{
    checkDefinition();
}


// Construct from components
Foam::faceCracker::faceCracker
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, Switch(dict.lookup("active"))),
    crackZoneID_
    (
        dict.lookup("faceZoneName"),
        mme.mesh().faceZones()
    ),
    crackPatchID_
    (
        dict.lookup("crackPatchName"),
        mme.mesh().boundaryMesh()
    ),
    coupledFacesToBreak_(),
    trigger_(false)
{
    checkDefinition();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceCracker::~faceCracker()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceCracker::setBreak
(
    const labelList& facesToBreak,
    const boolList& faceFlip,
    const labelList& coupledFacesToBreak
)
{
    if (trigger_)
    {
        FatalErrorIn
        (
            "void faceCracker::setBreak\n"
            "(\n"
            "    const labelList& facesToBreak,\n"
            "    const boolList& faceFlip,\n"
            "    const labelList& coupledFacesToBreak\n"
            ")"
        )   << "Setting faces to break before previous break was "
            << "not executed.  Probably an error in topo change handling."
            << abort(FatalError);
    }

    polyMesh& mesh = const_cast<polyMesh&>(topoChanger().mesh());

    // Check that all the faces in the face zone are internal
    if (debug)
    {
        // Check faces to break
        DynamicList<label> bouFacesInZone(facesToBreak.size());

        forAll (facesToBreak, faceI)
        {
            if (!mesh.isInternalFace(facesToBreak[faceI]))
            {
                bouFacesInZone.append(facesToBreak[faceI]);
            }
        }
    }

    // Put the faces into the face zone
    mesh.faceZones()[crackZoneID_.index()].resetAddressing
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


void Foam::faceCracker::setRefinement(polyTopoChange& ref) const
{
    if (trigger_)
    {
        //detachFaceCracker(ref);

        detachInternalFaces(ref);
        detachCoupledFaces(ref);

        // Reset the trigger
        trigger_ = false;
    }
}


void Foam::faceCracker::modifyMotionPoints(pointField& motionPoints) const
{}


void Foam::faceCracker::updateMesh(const mapPolyMesh&)
{
    // Mesh has changed topologically.  Update local topological data
    const polyMesh& mesh = topoChanger().mesh();

    crackZoneID_.update(mesh.faceZones());
    crackPatchID_.update(mesh.boundaryMesh());
}


void Foam::faceCracker::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name()<< nl
        << crackZoneID_.name() << nl
        << crackPatchID_.name() << endl;
}


void Foam::faceCracker::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << "    faceZoneName " << crackZoneID_.name()
        << token::END_STATEMENT << nl
        << "    crackPatchName " << crackPatchID_.name()
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
