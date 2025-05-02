/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "burstFvMeshTopoChanger.H"
#include "polyTopoChange.H"
#include "polyTopoChangeMap.H"
#include "volFields.H"
#include "PatchTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(burst, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, burst, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::burst::burst
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh)
{

    PtrList<entry> patches(dict.lookup("patches"));
    patchData_.setSize(patches.size());
    forAll(patches, i)
    {
        const dictionary& pDict = patches[i].dict();
        if
        (
            pDict.found("intactMasterPatch")
         && pDict.found("intactSlavePatch")
         && pDict.found("burstMasterPatch")
         && pDict.found("burstSlavePatch")
        )
        {
            patchData_.set
            (
                i,
                new burstData
                (
                    pDict.lookup<word>("intactMasterPatch"),
                    pDict.lookup<word>("intactSlavePatch"),
                    pDict.lookup<word>("burstMasterPatch"),
                    pDict.lookup<word>("burstSlavePatch"),
                    pDict
                )
            );
        }
        else if (pDict.found("intactPatch") && pDict.found("burstPatch"))
        {
            patchData_.set
            (
                i,
                new burstData
                (
                    pDict.lookup<word>("intactPatch"),
                    pDict.lookup<word>("burstPatch"),
                    pDict
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Couple burst patches require" << nl
                << "    intactMasterPatch" << nl
                << "    intactSlavePatch" << nl
                << "    burstMasterPatch" << nl
                << "    burstSlavePatch" << nl
                << nl
                << "Uncoupled burst patches require" << nl
                << "    intactPatch" << nl
                << "    burstPatch" << nl
                << endl
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::burst::~burst()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::burst::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
   DebugInfo<< "burst : Checking for topology changes..."
        << endl;

    DynamicList<label> masterFacesToChange;
    DynamicList<label> slaveFacesToChange;
    autoPtr<polyTopoChange> meshModPtr;
    forAll(patchData_, i)
    {
        masterFacesToChange.clear();
        slaveFacesToChange.clear();

        const burstData& data = patchData_[i];

        label burstMasterPatchID = - 1;
        label burstSlavePatchID = -1;

        label nMasterRemoved = 0;
        label nSlaveRemoved = 0;
        if (data.coupled())
        {
            burstMasterPatchID =
                mesh().boundaryMesh()[data.burstMasterPatch()].index();
            burstSlavePatchID =
                mesh().boundaryMesh()[data.burstSlavePatch()].index();

            const polyPatch& intactMasterPp =
                mesh().boundaryMesh()[data.intactMasterPatch()];
            const label intactMasterPatchID = intactMasterPp.index();

            const polyPatch& intactSlavePp =
                mesh().boundaryMesh()[data.intactSlavePatch()];
            const label intactSlavePatchID = intactSlavePp.index();

            data.burst().facesToChange
            (
                mesh().boundary()[intactMasterPatchID],
                mesh().boundary()[intactSlavePatchID],
                masterFacesToChange,
                slaveFacesToChange
            );

            nMasterRemoved = masterFacesToChange.size();
            nSlaveRemoved = slaveFacesToChange.size();
        }
        else
        {
            burstMasterPatchID =
                mesh().boundaryMesh()[data.burstPatch()].index();

            const polyPatch& intactPp =
                mesh().boundaryMesh()[data.intactPatch()];
            const label intactPatchID = intactPp.index();

            data.burst().facesToChange
            (
                mesh().boundary()[intactPatchID],
                masterFacesToChange
            );
            nMasterRemoved = masterFacesToChange.size();
        }
        reduce(nMasterRemoved, sumOp<label>());
        reduce(nSlaveRemoved, sumOp<label>());

        if (nMasterRemoved || nSlaveRemoved)
        {
            if (nMasterRemoved)
            {
                Info<< "Changing " << nMasterRemoved << " faces from "
                    << data.intactMasterPatch() << " to "
                    << data.burstMasterPatch() << endl;
            }
            if (nSlaveRemoved)
            {
                Info<< "Changing " << nSlaveRemoved << " faces from "
                    << data.intactSlavePatch() << " to "
                    << data.burstSlavePatch() << endl;
            }

            if (!meshModPtr.valid())
            {
                meshModPtr.set(new polyTopoChange(mesh()));
            }
            polyTopoChange& meshMod = meshModPtr();

            const faceList& faces = mesh().faces();
            const labelList& owner = mesh().faceOwner();
            forAll(masterFacesToChange, fi)
            {
                const label facei = masterFacesToChange[fi];
                meshMod.modifyFace
                (
                    faces[facei],
                    facei,
                    owner[facei],
                    -1,
                    false,
                    burstMasterPatchID
                );
            }
            forAll(slaveFacesToChange, fi)
            {
                const label facei = slaveFacesToChange[fi];
                meshMod.modifyFace
                (
                    faces[facei],
                    facei,
                    owner[facei],
                    -1,
                    false,
                    burstSlavePatchID
                );
            }
        }
    }

    if (meshModPtr.valid())
    {
        // Do any topology changes
        autoPtr<polyTopoChangeMap> map = meshModPtr->changeMesh(mesh());
        mesh().topoChange(map);
        return true;
    }

    return false;
}


void Foam::fvMeshTopoChangers::burst::topoChange(const polyTopoChangeMap& map)
{
    forAll(patchData_, i)
    {
        patchData_[i].burst().needUpdate();
    }
}


void Foam::fvMeshTopoChangers::burst::mapMesh(const polyMeshMap& map)
{
    forAll(patchData_, i)
    {
        patchData_[i].burst().needUpdate();
    }
}


void Foam::fvMeshTopoChangers::burst::distribute
(
    const polyDistributionMap& map
)
{
    forAll(patchData_, i)
    {
        patchData_[i].burst().needUpdate();
    }
}


// ************************************************************************* //
