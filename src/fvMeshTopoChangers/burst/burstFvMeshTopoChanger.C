/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
#include "forwardFieldMapper.H"
#include "fvMeshTools.H"
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

void Foam::fvMeshTopoChangers::burst::addProcessorCyclicPatches
(
    const List<bool>& hasOldSize
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    DynamicList<polyPatch*, 1, 1, 1> newPatches(mesh().boundaryMesh().size());

    label newPatchi = mesh().boundaryMesh().size();

    wordList addedPatches;

    // Add the processor cyclic patches
    forAll(patchData_, couplei)
    {
        const Pair<word> origPatchNames
        (
            patchData_[couplei].burstMasterPatch(),
            patchData_[couplei].burstSlavePatch()
        );
        Pair<word> ncPatchNames
        (
            nonConformalCyclicPolyPatch::typeName + "_on_" + origPatchNames[0],
            nonConformalCyclicPolyPatch::typeName + "_on_" + origPatchNames[1]
        );

        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        const polyPatch& patch1 = patches[origPatchNames.first()];
        const polyPatch& patch2 = patches[origPatchNames.second()];

        if
        (
            returnReduce
            (
                hasOldSize[couplei]

             && !patch1.size()
             && !patch2.size(),
                orOp<bool>()
            )
        )
        {
            continue;
        }

        boolList procHasPatch1(Pstream::nProcs(), false);
        procHasPatch1[Pstream::myProcNo()] = !patch1.empty();
        Pstream::gatherList(procHasPatch1);
        Pstream::scatterList(procHasPatch1);

        boolList procHasPatch2(Pstream::nProcs(), false);
        procHasPatch2[Pstream::myProcNo()] = !patch2.empty();
        Pstream::gatherList(procHasPatch2);
        Pstream::scatterList(procHasPatch2);

        // Multiple cyclic interfaces must be ordered in a specific way for
        // processor communication to function correctly.
        //
        // A communication that is sent from the cyclic owner is received
        // on the cyclic neighbour and vice versa. Therefore, in a coupled
        // pair of processors if one sends the owner first the other must
        // receive the neighbour first.
        //
        // We ensure the above by ordering the patches so that for the
        // lower indexed processor the owner interface comes first, and for
        // the higher indexed processor the neighbour comes first.

        auto appendProcPatches = [&](const bool owner, const bool first)
        {
            const boolList& procHasPatchA =
                owner ? procHasPatch1 : procHasPatch2;
            const boolList& procHasPatchB =
                owner ? procHasPatch2 : procHasPatch1;

            if (procHasPatchA[Pstream::myProcNo()])
            {
                forAll(procHasPatchB, proci)
                {
                    if
                    (
                        (
                            (first && proci > Pstream::myProcNo())
                         || (!first && proci < Pstream::myProcNo())
                        )
                     && procHasPatchB[proci]
                    )
                    {
                        autoPtr<nonConformalProcessorCyclicPolyPatch> ncpcpp
                        (
                            new nonConformalProcessorCyclicPolyPatch
                            (
                                0,
                                mesh().nFaces(),
                                patches.size(),
                                patches,
                                Pstream::myProcNo(),
                                proci,
                                ncPatchNames[!owner],
                                origPatchNames[!owner]
                            )
                        );

                        if (patches.findIndex(ncpcpp->name()) < 0)
                        {
                            addedPatches.append(ncpcpp->name());
                            newPatches(newPatchi++) = ncpcpp.ptr();
                        }
                    }
                }
            }
        };

        appendProcPatches(true, true);
        appendProcPatches(false, true);
        appendProcPatches(false, false);
        appendProcPatches(true, false);
    }

    if (returnReduce(newPatchi != mesh().boundaryMesh().size(), orOp<bool>()))
    {
        DebugInfo
            << "Adding " << nonConformalProcessorCyclicPolyPatch::typeName
            << " patches:" << nl
            << addedPatches << endl;

        const polyBoundaryMesh& patches = mesh().boundaryMesh();
        forAll(patches, patchi)
        {
            newPatches(patchi) = patches[patchi].clone(patches).ptr();
        }

        forAll(newPatches, newPatchi)
        {
            fvMeshTools::addPatch
            (
                mesh(),
                *newPatches[newPatchi],
                dictionary(),
                calculatedFvPatchField<scalar>::typeName,
                false
            );

            // Delete pointers
            deleteDemandDrivenData(newPatches[newPatchi]);
        }
    }
}


bool Foam::fvMeshTopoChangers::burst::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
   DebugInfo<< "burst : Checking for topology changes..."
        << endl;

    labelList boundaryMap(mesh().boundary().size(), -1);
    labelList rboundaryMap(mesh().boundary().size(), -1);
    List<Map<label>> bfaceMap(mesh().boundary().size());
    DynamicList<label> masterFacesToChange;
    DynamicList<label> slaveFacesToChange;
    autoPtr<polyTopoChange> meshModPtr;
    List<bool> hasOldSize(patchData_.size(), false);
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

            hasOldSize[i] = intactMasterPp.size() && intactSlavePp.size();

            data.burst().facesToChange
            (
                mesh().boundary()[intactMasterPatchID],
                mesh().boundary()[intactSlavePatchID],
                masterFacesToChange,
                slaveFacesToChange
            );

            nMasterRemoved = masterFacesToChange.size();
            nSlaveRemoved = slaveFacesToChange.size();

            if (nMasterRemoved)
            {
                boundaryMap[intactMasterPatchID] = burstMasterPatchID;
                rboundaryMap[burstMasterPatchID] = intactMasterPatchID;

                Map<label>& faceMap = bfaceMap[burstMasterPatchID];
                forAll(masterFacesToChange, fi)
                {
                    faceMap.insert(masterFacesToChange[fi], -1);
                }
            }
            if (nSlaveRemoved)
            {
                boundaryMap[intactSlavePatchID] = burstSlavePatchID;
                rboundaryMap[burstSlavePatchID] = intactSlavePatchID;

                Map<label>& faceMap = bfaceMap[burstSlavePatchID];
                forAll(slaveFacesToChange, fi)
                {
                    faceMap.insert(slaveFacesToChange[fi], -1);
                }
            }
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

            if (nMasterRemoved)
            {
                boundaryMap[intactPatchID] = burstMasterPatchID;
                rboundaryMap[burstMasterPatchID] = intactPatchID;

                Map<label>& faceMap = bfaceMap[burstMasterPatchID];
                forAll(masterFacesToChange, fi)
                {
                    faceMap.insert(masterFacesToChange[fi], -1);
                }
            }
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
        // Map all the volFields in the objectRegistry
        #define storeBoundariesType(Type, Mesh)                 \
            HashPtrTable<typename Mesh##Field<Type>::Boundary>  \
                bfields##Mesh##Type;                            \
            store##Mesh##Boundaries<Type>                       \
            (                                                   \
                boundaryMap,                                    \
                bfields##Mesh##Type                             \
            );
        FOR_ALL_FIELD_TYPES(storeBoundariesType, Vol);
        FOR_ALL_FIELD_TYPES(storeBoundariesType, Surface);

        // Do any topology changes
        autoPtr<polyTopoChangeMap> map = meshModPtr->changeMesh(mesh());
        mesh().topoChange(map);

        // Mapping from old patch to new patch using faceMap
        // New patch is created with out any source face so all values are Zero
        // Used copied patch from pre-mesh update to map unmapped faces
        // fvPatchFields are cloned so all data should be handled
        // Only new faces should be included so exisiting faces are left alone
        const labelList& faceMap = map().reverseFaceMap();
        const labelList& oldPatchStarts = map().oldPatchStarts();
        List<labelList> addressing(mesh().boundary().size());
        forAll(bfaceMap, patchi)
        {
            Map<label>& fMap = bfaceMap[patchi];
            const polyPatch& patch = mesh().boundaryMesh()[patchi];
            if (fMap.size())
            {
                const label start = patch.start();
                const label oldStart = oldPatchStarts[rboundaryMap[patchi]];

                labelList& addr = addressing[patchi];
                addr.setSize(patch.size(), -1);

                forAllIter(Map<label>, fMap, iter)
                {
                    const label oldLocalFacei = iter.key() - oldStart;
                    const label newLocalFacei = faceMap[iter.key()] - start;
                    addr[newLocalFacei] = oldLocalFacei;
                }
            }
        }


        bool needUpdate = false;
        forAll(hasOldSize, i)
        {
            if (!hasOldSize[i])
            {
                needUpdate = true;
            }
        }
        if (returnReduce(needUpdate, orOp<bool>()))
        {
            addProcessorCyclicPatches(hasOldSize);
        }

        // Map all the volFields in the objectRegistry using the mapping created
        // from the old local indices to the new local indices
        #define mapBoundariesType(Type, Mesh)                   \
            map##Mesh##Boundaries<Type>                         \
            (                                                   \
                addressing,                                     \
                bfields##Mesh##Type                             \
            );
        FOR_ALL_FIELD_TYPES(mapBoundariesType, Vol);
        FOR_ALL_FIELD_TYPES(mapBoundariesType, Surface);
        #undef mapBoundariesType

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
