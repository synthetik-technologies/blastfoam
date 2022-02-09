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
#include "fvMeshMapper.H"
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


wordList addEmptyPatches(fvMesh& mesh, const PtrList<entry>& entries)
{
    // Find new patches and add them to the mesh

    wordList addedPatches;
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

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
        &&  pbm.findIndex(patchName) < 0
        )
        {
            patchesToAdd.insert(patchName, patchDict);
            addedPatches.append(patchName);
        }
    }

    if (!patchesToAdd.size())
    {
        return addedPatches;
    }
    Info<< "Selected " << patchesToAdd.size() << " patches to add" << endl;

    // Add new patches with zero size to the mesh
    label curPatchIndex = pbm.size();
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
                curPatchIndex++,
                pbm
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

    // Make sure patches and zoneFaces are synchronised across couples
//     pbm.checkParallelSync(true);
//     mesh.faceZones().checkParallelSync(true);

    return addedPatches;
}


void updateTopoSets
(
    topoSetList& topoSets,
    PtrList<backupTopoSetSource>& regions,
    labelListList& selectedRegionCells,
    labelListList& selectedRegionFaces,
    boolListList& selectedRegionFlipMaps,
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
        boolList& selectedFlipMaps = selectedRegionFlipMaps[regionI];
        labelList& selectedPoints = selectedRegionPoints[regionI];

        selectedCells.clear();
        selectedFaces.clear();
        selectedFlipMaps.clear();
        selectedPoints.clear();
        regions[regionI].allowBackup(allowBackup);
        regions[regionI].createSets
        (
            selectedCells,
            selectedFaces,
            selectedFlipMaps,
            selectedPoints
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
            selectedFlipMaps,
            selectedPoints
        );
        Info<< endl;
    }
    topoSets.transferZones(false);
}


// Filter out the empty patches.
void filterPatches(fvMesh& mesh, const HashSet<word>& addedPatchNames)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
             || addedPatchNames.found(pp.name())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


void updateFields
(
    fvMesh& mesh,
    const mapPolyMesh& map,
    const HashSet<word>& addedPatches,
    const PtrList<entry>& entries
)
{
    // Correct boundary faces mapped-out-of-nothing.
    // This is just a hack to correct the value field.
    {
        fvMeshMapper mapper(mesh, map);
        bool hasWarned = false;

        forAllConstIter(HashSet<word>, addedPatches, iter)
        {
            label patchi = mesh.boundaryMesh().findPatchID(iter.key());

            const fvPatchMapper& pm = mapper.boundaryMap()[patchi];

            if (pm.sizeBeforeMapping() == 0)
            {
                if (!hasWarned)
                {
                    hasWarned = true;
                    WarningInFunction
                        << "Setting field on boundary faces to zero." << endl
                        << "You might have to edit these fields." << endl;
                }

                fvMeshTools::zeroPatchFields(mesh, patchi);
            }
        }
    }


    // Pass 2: change patchFields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        forAll(entries, entryI)
        {
            const dictionary& dict = entries[entryI].dict();
            if (dict.found("patches"))
            {
                const dictionary& patchSources = dict.subDict("patches");

                forAllConstIter(dictionary, patchSources, iter)
                {
                    const word patchName(iter().dict()["name"]);
                    label patchi = pbm.findPatchID(patchName);

                    if (iter().dict().found("patchFields"))
                    {
                        const dictionary& patchFieldsDict =
                            iter().dict().subDict
                            (
                                "patchFields"
                            );

                        fvMeshTools::setPatchFields
                        (
                            mesh,
                            patchi,
                            patchFieldsDict
                        );
                    }
                }
            }
            else
            {
                const dictionary& patchSource = dict.subDict("patchPairs");

                Switch sameGroup
                (
                    patchSource.lookupOrDefault("sameGroup", true)
                );

                const word& groupName = entries[entryI].name();

                if (patchSource.found("patchFields"))
                {
                    dictionary patchFieldsDict = patchSource.subDict
                    (
                        "patchFields"
                    );

                    if (sameGroup)
                    {
                        // Add coupleGroup to all entries
                        forAllIter(dictionary, patchFieldsDict, iter)
                        {
                            if (iter().isDict())
                            {
                                dictionary& dict = iter().dict();
                                dict.set("coupleGroup", groupName);
                            }
                        }

                        const labelList& patchIDs =
                            pbm.groupPatchIDs()[groupName];

                        forAll(patchIDs, i)
                        {
                            fvMeshTools::setPatchFields
                            (
                                mesh,
                                patchIDs[i],
                                patchFieldsDict
                            );
                        }
                    }
                    else
                    {
                        const word masterPatchName(groupName + "_master");
                        const word slavePatchName(groupName + "_slave");

                        label patchiMaster = pbm.findPatchID(masterPatchName);
                        label patchiSlave = pbm.findPatchID(slavePatchName);

                        fvMeshTools::setPatchFields
                        (
                            mesh,
                            patchiMaster,
                            patchFieldsDict
                        );

                        fvMeshTools::setPatchFields
                        (
                            mesh,
                            patchiSlave,
                            patchFieldsDict
                        );
                    }
                }
            }
        }
    }
}


// Synchronise points on both sides of coupled boundaries.
template<class CombineOp>
void syncPoints
(
    const polyMesh& mesh,
    pointField& points,
    const CombineOp& cop,
    const point& nullValue
)
{
    if (points.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << points.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if
            (
                isA<processorPolyPatch>(pp)
             && pp.nPoints() > 0
             && refCast<const processorPolyPatch>(pp).owner()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Get data per patchPoint in neighbouring point numbers.
                pointField patchInfo(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.nbrPoints();

                forAll(nbrPts, pointi)
                {
                    label nbrPointi = nbrPts[pointi];
                    if (nbrPointi >= 0 && nbrPointi < patchInfo.size())
                    {
                        patchInfo[nbrPointi] = points[meshPts[pointi]];
                    }
                }

                OPstream toNbr
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                );
                toNbr << patchInfo;
            }
        }


        // Receive and set.

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if
            (
                isA<processorPolyPatch>(pp)
             && pp.nPoints() > 0
             && !refCast<const processorPolyPatch>(pp).owner()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                pointField nbrPatchInfo(procPatch.nPoints());
                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }
                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

                if (procPatch.transform().transformsPosition())
                {
                    hasTransformation = true;
                    procPatch.transform().transformPosition
                    (
                        nbrPatchInfo,
                        nbrPatchInfo
                    );
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointi)
                {
                    label meshPointi = meshPts[pointi];
                    points[meshPointi] = nbrPatchInfo[pointi];
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            isA<cyclicPolyPatch>(pp)
         && refCast<const cyclicPolyPatch>(pp).owner()
        )
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();
            const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();
            const labelList& nbrMeshPts = nbrPatch.meshPoints();

            pointField patchPoints(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point0 = meshPts[e[0]];
                patchPoints[i] = points[point0];
            }

            if (cycPatch.transform().transformsPosition())
            {
                hasTransformation = true;
                cycPatch.transform().invTransformPosition
                (
                    patchPoints,
                    patchPoints
                );
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point1 = nbrMeshPts[e[1]];
                points[point1] = patchPoints[i];
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly necessary.
    // reduce(hasTransformation, orOp<bool>());

    // Synchronise multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }


        // Values on shared points.
        pointField sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = points[meshPointi];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            points[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
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
    boolListList savedFlipMaps(regions.size());
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
            savedFlipMaps,
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

    HashSet<word> addedPatches;

    // Find new patches and add them to the mesh
    if (setsToRemove.size())
    {
        Info<< "Removing cell sets" << endl;
        // Add new patches
        addedPatches.set(addEmptyPatches(mesh, setsToRemove));

        // Remove cells
        label nRemoved = 0;
        polyTopoChange meshMod(mesh);
        forAll(setsToRemove, seti)
        {
            const entry& patchEntry = setsToRemove[seti];
            const dictionary& dict = patchEntry.dict();
            const word patchName =
                dict.lookupOrDefault("patchName", patchEntry.keyword());

            Info<< "    Removing cells in " << patchEntry.keyword() << endl;
            nRemoved += meshTools::setRemoveCells
            (
                mesh,
                *topoSets.cellSets()[setsToRemove[seti].keyword()],
                patchName,
                meshMod,
                dict.lookupOrDefault("invert", false)
            );

        }
        Info<< "    Removed " << returnReduce(nRemoved, sumOp<label>())
            << " cells" << endl;
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());

        // Move mesh (since morphing might not do this)
//         if (map().hasMotionPoints())
//         {
//             mesh.movePoints(map().preMotionPoints());
//             mesh.moving(false);
//         }
    }

    // Find new patches and add them to the mesh
    if (patchesToAdd.size())
    {
        Info<< nl << "Adding patches" << endl;
        // Add new patches
        addedPatches.set(addEmptyPatches(mesh, patchesToAdd));

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

        // Move mesh (since morphing might not do this)
//         if (map().hasMotionPoints())
//         {
//             mesh.movePoints(map().preMotionPoints());
//             mesh.moving(false);
//         }
    }

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
        addedPatches.set(addEmptyPatches(mesh, bafflePatchesToAdd));

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

        // Move mesh (since morphing might not do this)
//         if (map().hasMotionPoints())
//         {
//             mesh.movePoints(map().preMotionPoints());
//             mesh.moving(false);
//         }
    }

    // Synchronise points.
//     if (!pointSync)
//     {
//         Info<< "Not synchronising points." << nl << endl;
//     }
//     else
//     {
//         Info<< "Synchronising points." << nl << endl;
//
//         pointField newPoints(mesh.points());
//
//         syncPoints
//         (
//             mesh,
//             newPoints,
//             minMagSqrEqOp<vector>(),
//             point(great, great, great)
//         );
//
//         scalarField diff(mag(newPoints-mesh.points()));
//         Info<< "Points changed by average:" << gAverage(diff)
//             << " max:" << gMax(diff) << nl << endl;
//
//         mesh.movePoints(newPoints);
//         mesh.moving(false);
//     }

//     autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
//     mesh.updateMesh(map());
//
//     // Move mesh (since morphing might not do this)
//     if (map().hasMotionPoints())
//     {
//         mesh.movePoints(map().preMotionPoints());
//     }


    if (overwrite)
    {
        mesh.setInstance(runTime.constant());
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
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh
            ),
            mesh.points()
        );
        points0.write();
        writeMesh = true;
    }

    // Remove any now zero-sized patches
    filterPatches(mesh, addedPatches);

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
