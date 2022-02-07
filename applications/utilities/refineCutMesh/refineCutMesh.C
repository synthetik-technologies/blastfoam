/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

Description
    Set values on a selected set of cells/patchfaces through a dictionary and
    refines using hexRef method.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "backupTopoSetSource.H"
#include "topoSetList.H"
#include "volFields.H"
#include "systemDict.H"
#include "PtrListDictionary.H"

#include "fvMeshHexRefiner.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "IOobjectList.H"
#include "dictionaryEntry.H"

#include "fvMeshTools.H"
#include "dynMeshTools.H"

using namespace Foam;


//- Read and add fields to the database
template<class Type, template<class> class Patch, class Mesh>
void readGeoFields(const fvMesh& mesh, const IOobjectList& objects)
{
    typedef GeometricField<Type, Patch, Mesh> FieldType;

    IOobjectList fields = objects.lookupClass(FieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (!mesh.foundObject<FieldType>(fieldIter()->name()))
        {
            IOobject fieldTargetIOobject
            (
                fieldIter()->name(),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldTargetIOobject.typeHeaderOk<FieldType>(true))
            {
                FieldType* fPtr
                (
                    new FieldType
                    (
                        fieldTargetIOobject,
                        mesh
                    )
                );
                fPtr->store(fPtr);
            }
        }
    }
}


//- Read and add fields to the database
template<class Type>
void readPointFields(const fvMesh& mesh, const IOobjectList& objects)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> FieldType;
    IOobjectList fields = objects.lookupClass(FieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (!mesh.foundObject<FieldType>(fieldIter()->name()))
        {
            IOobject fieldTargetIOobject
            (
                fieldIter()->name(),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldTargetIOobject.typeHeaderOk<FieldType>(true))
            {
                FieldType* fPtr
                (
                    new FieldType
                    (
                        fieldTargetIOobject,
                        pointMesh::New(mesh)
                    )
                );
                fPtr->store(fPtr);
            }
        }
    }
}


//- Read and add all fields to the database
void readAndAddAllFields(const fvMesh& mesh)
{
    // Get all fields present at the current time
    IOobjectList objects(mesh, mesh.time().timeName());

    readGeoFields<scalar, fvPatchField, volMesh>(mesh, objects);
    readGeoFields<vector, fvPatchField, volMesh>(mesh, objects);
    readGeoFields<symmTensor, fvPatchField, volMesh>(mesh, objects);
    readGeoFields<sphericalTensor, fvPatchField, volMesh>(mesh, objects);
    readGeoFields<tensor, fvPatchField, volMesh>(mesh, objects);

    readGeoFields<scalar, fvsPatchField, surfaceMesh>(mesh, objects);
    readGeoFields<vector, fvsPatchField, surfaceMesh>(mesh, objects);
    readGeoFields<symmTensor, fvsPatchField, surfaceMesh>(mesh, objects);
    readGeoFields<sphericalTensor, fvsPatchField, surfaceMesh>(mesh, objects);
    readGeoFields<tensor, fvsPatchField, surfaceMesh>(mesh, objects);

    readPointFields<scalar>(mesh, objects);
    readPointFields<vector>(mesh, objects);
    readPointFields<symmTensor>(mesh, objects);
    readPointFields<sphericalTensor>(mesh, objects);
    readPointFields<tensor>(mesh, objects);
}


void addEmptyPatches(fvMesh& mesh, const PtrList<entry>& entries)
{
    // Find new patches and add them to the mesh

    // Collect all new patches
    HashTable<dictionary> patchesToAdd;
    forAll(entries, seti)
    {
        const entry& patchEntry = entries[seti];
        const dictionary& patchDict(patchEntry.dict());
        const word patchName = patchDict.lookupOrDefault<word>
        (
            "patchName",
            patchEntry.keyword()
        );
        if
        (
            !patchesToAdd.found(patchName)
        &&  mesh.boundaryMesh().findIndex(patchName) < 0
        )
        {
            patchesToAdd.insert(patchName, patchDict);
        }
    }

    if (!patchesToAdd.size())
    {
        return;
    }
    Info<< "Selected " << patchesToAdd.size() << " patches to add" << endl;

    // Add new patches with zero size to the mesh
    label curPatchIndex = mesh.boundaryMesh().size();
    forAllConstIter(HashTable<dictionary>, patchesToAdd, iter)
    {
        dictionary dict(iter());
        dict.set("nFaces", 0);
        dict.set("startFace", 0);
        word groupName = word::null;
        if (dict.found("group"))
        {
            groupName = dict.lookup<word>("group");
            dict.remove("group");
        }
        autoPtr<polyPatch> ppPtr
        (
            polyPatch::New
            (
                iter.key(),
                dict,
                curPatchIndex,
                mesh.boundaryMesh()
            )
        );
        polyPatch& pp = ppPtr();
        if (!groupName.empty() && !pp.inGroup(groupName))
        {
            pp.inGroups().append(groupName);
        }

        Info<< "    Adding new patch " << iter.key() << endl;
        fvMeshTools::addPatch
        (
            mesh,
            pp,
            dictionary(),
            calculatedFvPatchField<scalar>::typeName,
            true
        );
    }
}


void updateTopoSets
(
    topoSetList& topoSets,
    const PtrList<backupTopoSetSource>& regions,
    labelListList& selectedRegionCells,
    labelListList& selectedRegionFaces,
    labelListList& selectedRegionPoints,
    const bool allowBackup = false
)
{
    Info<< "Setting field region values" << endl;
    forAll(regions, regionI)
    {
        const dictionary& regionDict =  regions[regionI].dict();

        labelList& selectedCells = selectedRegionCells[regionI];
        labelList& selectedFaces = selectedRegionFaces[regionI];
        labelList& selectedPoints = selectedRegionPoints[regionI];

        selectedCells.clear();
        selectedFaces.clear();
        selectedPoints.clear();
        regions[regionI].createSets
        (
            selectedCells,
            selectedFaces,
            selectedPoints,
            allowBackup
        );

        if
        (
            regions[regionI].isCell()
            && !returnReduce(selectedCells.size(), sumOp<label>())
        )
        {
            WarningInFunction
                << "No cells were selected for using " << nl
                << regionDict
                << "To expand searchable region add backup " << nl
                << "Region or expand backup region." << endl;
        }
        else if
        (
            regions[regionI].isFace()
            && !returnReduce(selectedFaces.size(), sumOp<label>())
        )
        {
            WarningInFunction
                << "No faces were selected for using " << nl
                << regionDict
                << "To expand searchable region add backup " << nl
                << "Region or expand backup region." << endl;
        }
        else if
        (
            regions[regionI].isPoint()
            && !returnReduce(selectedPoints.size(), sumOp<label>())
        )
        {
            WarningInFunction
                << "No points were selected for using " << nl
                << regionDict
                << "To expand searchable region add backup " << nl
                << "Region or expand backup region." << endl;
        }

        bool set =
            selectedCells.size()
         || selectedFaces.size()
         || selectedPoints.size();
        reduce(set, orOp<bool>());

        if (set)
        {
            Info<< "    Selected "
                << returnReduce(selectedCells.size(), sumOp<label>())
                << " cells, "
                << returnReduce(selectedFaces.size(), sumOp<label>())
                << " faces, "
                << returnReduce(selectedPoints.size(), sumOp<label>())
                << " points" << endl;
        }

        topoSets.update
        (
            regionDict,
            selectedCells,
            selectedFaces,
            selectedPoints
        );
        Info<< endl;
    }
    topoSets.transferZones(false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //- Add options
    timeSelector::addOptions(true, false);

    argList::addBoolOption
    (
        "noFields",
        "Do not update fields"
    );
    argList::addBoolOption
    (
        "debug",
        "Output partial updates and additional fields"
    );
    argList::addBoolOption
    (
        "forceHex8",
        "Force use of standard OpenFOAM hexRef8 refinement"
    );
    argList::addBoolOption
    (
        "noRefine",
        "Do not refine"
    );
    argList::addBoolOption
    (
        "overwrite",
        "Write the mesh to constant"
    );
    argList::addBoolOption
    (
        "noHistory",
        "Do not write the history"
    );

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Select time
    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const dictionary refineCutDict(systemDict("refineCutDict", args, mesh));

    // Do not load in fields
    bool noFields(args.optionFound("noFields"));

    // Do not write fields
    // Usefull if a refined mesh is needed before mesh manipulation
    bool noWrite(args.optionFound("noWrite") || noFields);

    bool noRefine(args.optionFound("noRefine"));
    bool overwrite(args.optionFound("overwrite"));
    bool noHistory(args.optionFound("noHistory"));
    bool debug(args.optionFound("debug"));

    //- Is the mesh balanced
    autoPtr<fvMeshRefiner> refiner;
    if (!noRefine)
    {
        dictionary refineDict
        (
            refineCutDict.optionalSubDict
            (
                "refinerCoeffs"
            )
        );
        if (args.optionFound("forceHex8"))
        {
            refineDict.set("forceHex8", true);
            refiner.set(new fvMeshHexRefiner(mesh, refineDict, true));
        }
        else
        {
            refiner = fvMeshRefiner::New(mesh, refineDict, true);
        }
        refiner->setForce(true);
    }
    bool refine = refiner.valid();

    volScalarField error
    (
        IOobject
        (
            "error",
            runTime.timeName(),
            mesh
        ),
        mesh,
        0.0
    );

    // Read in all fields to allow resizing
    if (!noFields)
    {
        readAndAddAllFields(mesh);
    }

    //- List of sources (and backups if present)
    // Stored to reduce the number of reads
    PtrList<backupTopoSetSource> regions
    (
        refineCutDict.lookup("regions"),
        backupTopoSetSource::iNew(mesh)
    );
    topoSetList topoSets(mesh);

    labelList levels(regions.size(), 0);
    forAll(regions, regionI)
    {
        levels[regionI] =
            regions[regionI].dict().lookupOrDefault("level", 0);
    }

    label maxLevel =
    (
        levels.size() && !refineCutDict.found("maxRefinement")
      ? max(levels)
      : refineCutDict.lookupOrDefault<label>("maxRefinement", 0)
    );


    // Maximum number of iterations
    label iter = 0;
    label maxIter = 1;
    if (refine)
    {
        maxIter =
            max
            (
                1,
                refineCutDict.lookupOrDefault
                (
                    "maxIter",
                    max
                    (
                        2*maxLevel,
                        gMax(refiner->cellLevel())*2
                    )
                )
            );
    }

    // Flag for final iteration
    bool end = false;

    // Flag to initiate end
    bool prepareToStop = (maxIter == 1);

    // List of saved cells (per cell set)
    labelListList savedCells(regions.size());
    labelListList savedFaces(regions.size());
    labelListList savedPoints(regions.size());

    while(!end)
    {
        if (debug)
        {
            runTime++;
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        if (maxIter <= iter)
        {
            prepareToStop = true;
        }

        error = -1.0;

        // Check if this is the final iteration so correct cell shapes are used
        if (prepareToStop)
        {
            end = true;
        }

        updateTopoSets
        (
            topoSets,
            regions,
            savedCells,
            savedFaces,
            savedPoints,
            !end
        );

        // Update error and mesh if not the final iteration
        if (refine)
        {
            labelList maxCellLevel(mesh.nCells(), -1);
            forAll(regions, regionI)
            {
                // Set specified cells to be refined
                labelList refineCells;
                labelList refineFaces;
                labelList refinePoints;
                if
                (
                    regions[regionI].dict().lookupOrDefault
                    (
                        "refineInternal",
                        false
                    )
                )
                {
                    refineCells = savedCells[regionI];
                }
                else if
                (
                    regions[regionI].dict().lookupOrDefault
                    (
                        "refineInterface",
                        false
                    )
                )
                {
                    refineCells =
                        topoSets.extractInterfaceCells(savedCells[regionI]);
                }
                if
                (
                    regions[regionI].dict().lookupOrDefault
                    (
                        "refineFaces",
                        false
                    )
                )
                {
                    refineFaces = topoSets.extractSelectedFaces
                    (
                        regions[regionI].dict(),
                        savedFaces[regionI],
                        true
                    );
                }
                if
                (
                    regions[regionI].dict().lookupOrDefault
                    (
                        "refinePoints",
                        false
                    )
                )
                {
                    refinePoints = topoSets.extractSelectedPoints
                    (
                        regions[regionI].dict(),
                        savedPoints[regionI],
                        true
                    );
                }

                // Do not use the max level, use current
                // Order is important in the definitions of regions
                Switch overwriteLevel =
                    regions[regionI].dict().lookupOrDefault<Switch>
                    (
                        "overwriteLevel",
                        false
                    );

                // Set actual max cell level
                const labelList& owner = mesh.faceOwner();
                const labelList& neighbour = mesh.faceNeighbour();
                const labelListList& pointCells = mesh.pointCells();
                if (overwriteLevel)
                {
                    forAll(refineCells, celli)
                    {
                        maxCellLevel[refineCells[celli]] = levels[regionI];
                    }

                    forAll(refineFaces, fi)
                    {
                        const label facei = refineFaces[fi];
                        maxCellLevel[owner[facei]] = levels[regionI];
                        if (facei < mesh.nInternalFaces())
                        {
                            maxCellLevel[neighbour[facei]] = levels[regionI];
                        }
                    }
                    forAll(refinePoints, pi)
                    {
                        const label pointi = refinePoints[pi];
                        const labelList& pc = pointCells[pointi];
                        forAll(pc, ci)
                        {
                            maxCellLevel[pc[ci]] = levels[regionI];
                        }
                    }
                }
                else
                {
                    forAll(refineCells, celli)
                    {
                        maxCellLevel[refineCells[celli]] =
                            max
                            (
                                maxCellLevel[refineCells[celli]],
                                levels[regionI]
                            );
                    }
                    forAll(refineFaces, fi)
                    {
                        const label facei = refineFaces[fi];
                        maxCellLevel[owner[facei]] =
                            max
                            (
                                maxCellLevel[owner[facei]],
                                levels[regionI]
                            );
                        if (facei < mesh.nInternalFaces())
                        {
                            maxCellLevel[neighbour[facei]] =
                                max
                                (
                                    maxCellLevel[neighbour[facei]],
                                    levels[regionI]
                                );
                        }
                    }
                    forAll(refinePoints, pi)
                    {
                        const label pointi = refinePoints[pi];
                        const labelList& pc = pointCells[pointi];
                        forAll(pc, ci)
                        {
                            maxCellLevel[pc[ci]] =
                                max
                                (
                                    maxCellLevel[pc[ci]],
                                    levels[regionI]
                                );
                        }
                    }
                }

                forAll(refineCells, celli)
                {
                    error[refineCells[celli]] = 1.0;
                }
                forAll(refineFaces, fi)
                {
                    const label facei = refineFaces[fi];
                    error[owner[facei]] = 1.0;

                    if (facei < mesh.nInternalFaces())
                    {
                        error[neighbour[facei]] = 1.0;
                    }
                }
                forAll(refinePoints, pi)
                {
                    const label pointi = refinePoints[pi];
                    const labelList& pc = pointCells[pointi];
                    forAll(pc, ci)
                    {
                        error[pc[ci]] = 1.0;
                    }
                }

                // Extend refinement by nBufferLayers
                for
                (
                    label i = 0;
                    i < refiner->nRefinementBufferLayers() + 1;
                    i++
                )
                {
                    fvMeshRefiner::extendMaxCellLevel
                    (
                        mesh,
                        savedCells[regionI],
                        maxCellLevel,
                        levels[regionI]
                    );
                }
            }

            labelList maxRefinement(mesh.nCells(), maxLevel);

            // Set the maxCell level
            forAll(maxCellLevel, celli)
            {
                if (maxCellLevel[celli] < 0)
                {
                    maxCellLevel[celli] = maxRefinement[celli];
                }
            }

            // Mark cells greater than the max cell level for unrefinment
            const labelList& cellLevel = refiner->cellLevel();
            forAll(error, celli)
            {
                if (cellLevel[celli] == maxCellLevel[celli])
                {
                    error[celli] = 0.0;
                }
                else if (cellLevel[celli] > maxCellLevel[celli])
                {
                    error[celli] = -1.0;
                }
            }

            // Write fields and mesh if using debug
            if (debug)
            {
                mesh.setInstance(runTime.timeName());
                bool writeOk = (mesh.write() && refiner->write());
                volScalarField scalarMaxCellLevel
                (
                    volScalarField::New
                    (
                        "maxCellLevel",
                        mesh,
                        dimensionedScalar(dimless, 0),
                        extrapolatedCalculatedFvPatchField<scalar>::typeName
                    )
                );

                forAll(cellLevel, celli)
                {
                    scalarMaxCellLevel[celli] = maxCellLevel[celli];
                }
                scalarMaxCellLevel.correctBoundaryConditions();
                writeOk =
                    writeOk
                 && scalarMaxCellLevel.write();

                Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
            }

            // Update mesh (return if mesh changes)
            if (!end)
            {
                prepareToStop = !refiner->refine(error, maxCellLevel);
            }
        }
        iter++;
    }

    if (noHistory)
    {
        refiner.clear();
    }


    // Remove cells in the selected regions
    PtrList<entry> setsToRemove(refineCutDict.lookup("remove"));

    // Remove cells in the selected regions
    PtrList<entry> patchesToAdd(refineCutDict.lookup("patches"));

    // Add baffles
    PtrList<entry> bafflesToAdd(refineCutDict.lookup("baffles"));
    PtrListDictionary<entry> bafflePatchesToAdd(bafflesToAdd.size()*2);
    List<Pair<word>> bafflePatches(bafflesToAdd.size());


    // Find new patches and add them to the mesh
    if (setsToRemove.size())
    {
        Info<< "Removing cell sets" << endl;
        // Add new patches
        addEmptyPatches(mesh, setsToRemove);

        // Remove cells
        polyTopoChange meshMod(mesh);
        forAll(setsToRemove, seti)
        {
            const entry& patchEntry = setsToRemove[seti];
            const dictionary& dict = patchEntry.dict();
            const word patchName =
                dict.lookupOrDefault("patchName", patchEntry.keyword());

            Info<< "    Removing cells in " << patchEntry.keyword() << endl;
            meshTools::setRemoveCells
            (
                mesh,
                *topoSets.cellSets()[setsToRemove[seti].keyword()],
                patchName,
                meshMod,
                dict.lookupOrDefault("invert", false)
            );
        }
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
    }

//     updateTopoSets
//     (
//         topoSets,
//         regions,
//         savedCells,
//         savedFaces,
//         savedPoints
//     );

    // Find new patches and add them to the mesh
    if (patchesToAdd.size())
    {
        Info<< nl << "Adding patches" << endl;
        // Add new patches
        addEmptyPatches(mesh, patchesToAdd);

        // Modify faces
        polyTopoChange meshMod(mesh);
        forAll(patchesToAdd, seti)
        {
            const entry& patchEntry = patchesToAdd[seti];
            const dictionary& patchDict = patchEntry.dict();
            const word patchName = patchEntry.keyword();
            const word setName = patchDict.lookup<word>("set");
            const topoSet& set = *topoSets.faceSets()[setName];
            PackedBoolList modifiedFace(mesh.nFaces());

            faceZone fZone
            (
                setName,
                set.toc(),
                boolList(set.size(), false),
                mesh.faceZones().findIndex(setName),
                mesh.faceZones()
            );

            meshTools::createPatchFaces
            (
                false,
                mesh,
                fZone,
                {mesh.boundaryMesh().findIndex(patchName)},
                meshMod,
                modifiedFace
            );
        }
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
    }

//     updateTopoSets
//     (
//         topoSets,
//         regions,
//         savedCells,
//         savedFaces,
//         savedPoints
//     );

    // Find new patches and add them to the mesh
    if (bafflesToAdd.size())
    {
        Info<< "Adding Baffles" << endl;
        forAll(bafflesToAdd, bafflei)
        {
            const entry& baffleEntry = bafflesToAdd[bafflei];
            dictionary baffleDict(baffleEntry.dict());
            const word groupName =
                baffleDict.lookupOrDefault("patchName", baffleEntry.keyword());
            const word masterName
            (
                baffleDict.lookupOrDefault<word>
                (
                    "masterPatch",
                    groupName + "_master"
                )
            );
            const word slaveName
            (
                baffleDict.lookupOrDefault<word>
                (
                    "slavePatch",
                    groupName + "_slave"
                )
            );
            if (baffleDict.found("patchName"))
            {
                baffleDict.remove("patchName");
            }
            baffleDict.remove("set");
            dictionary masterDict
            (
                baffleDict.isDict(masterName)
              ? baffleDict.subDict(masterName)
              : baffleDict
            );
            dictionary slaveDict
            (
                baffleDict.isDict(slaveName)
              ? baffleDict.subDict(slaveName)
              : baffleDict
            );

            Switch sameGroup
            (
                baffleDict.lookupOrDefault("sameGroup", true)
            );
            word masterGroupName = groupName;
            word slaveGroupName = groupName;
            if (!sameGroup)
            {
                masterGroupName = groupName + "_master";
                slaveGroupName = groupName + "_slave";
            }
            masterDict.set("coupleGroup", masterGroupName);
            masterDict.set("group", masterGroupName);
            slaveDict.set("coupleGroup", slaveGroupName);
            slaveDict.set("group", slaveGroupName);


            bafflePatchesToAdd.set
            (
                2*bafflei,
                masterName,
                new dictionaryEntry
                (
                    masterName,
                    baffleEntry.dict(),
                    masterDict
                )
            );
            bafflePatchesToAdd.set
            (
                2*bafflei + 1,
                slaveName,
                new dictionaryEntry
                (
                    slaveName,
                    baffleEntry.dict(),
                    slaveDict
                )
            );
            bafflePatches[bafflei] = Pair<word>(masterName, slaveName);
        }

        // Add new patches
        addEmptyPatches(mesh, bafflePatchesToAdd);

        // Modify faces
        polyTopoChange meshMod(mesh);
        forAll(bafflesToAdd, bafflei)
        {
            const entry& patchEntry = bafflesToAdd[bafflei];
            const dictionary& patchDict = patchEntry.dict();
            const word patchName = patchEntry.keyword();
            const word setName = patchDict.lookup<word>("set");
            PackedBoolList modifiedFace(mesh.nFaces());
            const label zoneID = mesh.faceZones().findIndex(setName);

            const faceZone& fZone(mesh.faceZones()[zoneID]);
            Info<<returnReduce(fZone.size(), sumOp<label>())<<endl;

            meshTools::createBaffleFaces
            (
                false,
                mesh,
                fZone,
                {mesh.boundaryMesh().findIndex(bafflePatches[bafflei][0])},
                {mesh.boundaryMesh().findIndex(bafflePatches[bafflei][1])},
                meshMod,
                modifiedFace
            );
        }
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
    }


    bool writeMesh = topoSets.writeSets();
    if (refine && !debug)
    {
        //- Write points0 field to time directory
        pointIOField points0
        (
            IOobject
            (
                "points0",
                overwrite ? runTime.constant() : runTime.timeName(),
                polyMesh::meshSubDir,
                mesh
            ),
            mesh.points()
        );
        points0.write();
        writeMesh = true;
    }


//     if (writeMesh)
    {
        mesh.write();
    }

    Info<< "\nEnd\n" << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
