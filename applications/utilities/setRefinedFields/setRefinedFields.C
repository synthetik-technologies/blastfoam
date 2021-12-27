/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-10-2019  Jeff Heylmun:   Added refinement to setFields utility
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
#include "topoSetSource.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "volFields.H"
#include "systemDict.H"

#include "hexRef.H"
#include "hexRef3D.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"
#include "errorEstimator.H"

#include "IOobjectList.H"
#include "FieldSetType.H"

using namespace Foam;

#include "refineFunctions.H"


enum GeoType
{
    VOL,
    SURFACE,
    POINT,
    UNKNOWN_GEO
};
HashTable<GeoType> geoTypeNames
(
    {
        {"vol", VOL},
        {"surface", SURFACE},
        {"point", POINT}
    }
);
HashTable<word, GeoType, Hash<label>> geoEnumTypes
(
    {
        {VOL, "vol"},
        {SURFACE, "surface"},
        {POINT, "point"}
    }
);


enum PrimitiveType
{
    SCALAR,
    VECTOR,
    SPHERICALTENSOR,
    SYMMTENSOR,
    TENSOR,
    UNKNOWN_PRIM
};
HashTable<PrimitiveType> primitiveTypeNames
(
    {
        {"Scalar", SCALAR},
        {"Vector", VECTOR},
        {"SymmTensor", SYMMTENSOR},
        {"SphericalTensor", SPHERICALTENSOR},
        {"Tensor", TENSOR}
    }
);
HashTable<word, PrimitiveType, Hash<label>> primitiveEnumTypes
(
    {
        {SCALAR, "vol"},
        {VECTOR, "surface"},
        {SYMMTENSOR, "SymmTensor"},
        {SPHERICALTENSOR, "SphericalTensor"},
        {TENSOR, "Tensor"}
    }
);

GeoType getGeoType(const word& type)
{
    forAllConstIter(HashTable<GeoType>, geoTypeNames, iter)
    {
        if (label(type.find(iter.key())) >= 0)
        {
            return iter();
        }
    }
//     FatalErrorInFunction
//         << "Could not determine geometry type from " << type << endl
//         << abort(FatalError);
    return UNKNOWN_GEO;
}


PrimitiveType getPrimitiveType(const word& type)
{
    forAllConstIter(HashTable<PrimitiveType>, primitiveTypeNames, iter)
    {
        if (label(type.find(iter.key())) >= 0)
        {
            return iter();
        }
    }
//     FatalErrorInFunction
//         << "Could not determine primitive type from " << type << endl
//         << abort(FatalError);
    return UNKNOWN_PRIM;
}


//- Select faces based on selected cells
labelList cellsToFaces(const fvMesh& mesh, const labelList& cells)
{
    labelHashSet selectedFaces;
    forAll(cells, ci)
    {
        selectedFaces.insert(mesh.cells()[cells[ci]]);
    }
    return selectedFaces.toc();
}

//- Select Points based on selected cells
labelList cellsToPoints(const fvMesh& mesh, const labelList& cells)
{
    labelHashSet selectedPoints;
    forAll(cells, ci)
    {
        selectedPoints.insert(mesh.cellPoints()[cells[ci]]);
    }
    return selectedPoints.toc();
}

//- Select Cells based on selected faces
labelList facesToCells(const fvMesh& mesh, const labelList& faces)
{
    if (!faces.size())
    {
        return labelList();
    }
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    labelHashSet selectedFaces(faces);
    labelHashSet selectedCells;
    forAll(owner, facei)
    {
        if (selectedFaces.found(owner[facei]))
        {
            if (facei >= mesh.nInternalFaces())
            {
                selectedCells.insert(owner[facei]);
            }
            else if (selectedFaces.found(neighbour[facei]))
            {
                selectedCells.insert(owner[facei]);
            }
        }
    }
    return selectedCells.toc();
}

//- Select Points based on selected faces
labelList facesToPoints(const fvMesh& mesh, const labelList& faces)
{
    labelHashSet selectedPoints;
    forAll(faces, fi)
    {
        selectedPoints.insert(mesh.faces()[faces[fi]]);
    }
    return selectedPoints.toc();
}

//- Select cells based on selected points
labelList pointsToCells(const fvMesh& mesh, const labelList& points)
{
    if (!points.size())
    {
        return labelList();
    }
    labelHashSet selectedPoints(points);
    labelHashSet selectedCells;
    forAll(mesh.cells(), celli)
    {
        const cell& c = mesh.cells()[celli];
        const labelList cp = c.labels(mesh.faces());
        bool allFound = true;
        forAll(cp, pi)
        {
            if (!selectedPoints.found(cp[pi]))
            {
                allFound = false;
                break;
            }
        }
        if (allFound)
        {
            selectedCells.insert(celli);
        }
    }
    return selectedCells.toc();
}

//- Select faces based on selected points
labelList pointsToFaces(const fvMesh& mesh, const labelList& points)
{
    if (!points.size())
    {
        return labelList();
    }
    labelHashSet selectedPoints(points);
    labelHashSet selectedFaces;
    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];
        bool allFound = true;
        forAll(f, pi)
        {
            if (!selectedPoints.found(f[pi]))
            {
                allFound = false;
                break;
            }
        }
        if (allFound)
        {
            selectedFaces.insert(facei);
        }
    }
    return selectedFaces.toc();
}

class topoSetList
{

public:

    topoSetList()
    {}

    autoPtr<topoSetList> clone() const
    {
        return autoPtr<topoSetList>(new topoSetList());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const dictionary& dict_;
        const labelList selectedCells_;
        const labelList selectedFaces_;
        const labelList selectedPoints_;
        const bool write_;
        const bool force_;

        template<class Type, template<class> class Patch, class Mesh>
        void createTopoSetType
        (
            const word& fieldSetTypeDesc,
            const labelList& elms,
            Istream& is
        ) const
        {
            token t(is);
            if (!t.isWord())
            {
                return;
            }
            const word fieldName(t.wordToken());
            autoPtr<FieldSetType<Type, Patch, Mesh>> fieldSet
            (
                FieldSetType<Type, Patch, Mesh>::New
                (
                    fieldSetTypeDesc,
                    fieldName,
                    mesh_,
                    dict_,
                    elms,
                    is,
                    write_
                )
            );
            if (fieldSet->good())
            {
                Info<< "    Setting " << fieldName << endl;
            }
            else
            {
                WarningInFunction
                    << "Field " << fieldName << " not found" << endl;
            }
        }

        template<template<class> class Patch, class Mesh>
        void createTopoSet
        (
            const word& fieldSetType,
            const labelList& elms,
            Istream& is
        ) const
        {
            if (!elms.size() && !force_)
            {
                return;
            }

            PrimitiveType prim(getPrimitiveType(fieldSetType));
            switch (prim)
            {
                case SCALAR:
                    createTopoSetType<scalar, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case VECTOR:
                    createTopoSetType<vector, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case SYMMTENSOR:
                    createTopoSetType<symmTensor, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case SPHERICALTENSOR:
                    createTopoSetType<sphericalTensor, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                case TENSOR:
                    createTopoSetType<tensor, Patch, Mesh>
                    (
                        fieldSetType,
                        elms,
                        is
                    );
                    break;
                default:
                    break;
            }
        }


    public:

        iNew
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const bool force = true
        )
        :
            mesh_(mesh),
            dict_(dict),
            write_(false),
            force_(force)
        {}

        iNew
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelList& selectedCells,
            const labelList& selectedFaces,
            const labelList& selectedPoints,
            const bool write,
            const bool force = false
        )
        :
            mesh_(mesh),
            dict_(dict),
            selectedCells_(selectedCells),
            selectedFaces_(selectedFaces),
            selectedPoints_(selectedPoints),
            write_(write),
            force_(force)
        {}

        autoPtr<topoSetList> operator()(Istream& is) const
        {
            token t(is);
            if (!t.isWord())
            {
                return autoPtr<topoSetList>();
            }

            word fieldSetType(t.wordToken());
            GeoType geo(getGeoType(fieldSetType));

            switch (geo)
            {
                case VOL:
                    createTopoSet<fvPatchField, volMesh>
                    (
                        fieldSetType,
                        selectedCells_,
                        is
                    );
                    break;
                case SURFACE:
                    createTopoSet<fvsPatchField, surfaceMesh>
                    (
                        fieldSetType,
                        selectedFaces_,
                        is
                    );
                    break;
                case POINT:
                    createTopoSet<pointPatchField, pointMesh>
                    (
                        fieldSetType,
                        selectedPoints_,
                        is
                    );
                    break;
                default:
                    break;
            }

            return autoPtr<topoSetList>(new topoSetList());
        }
    };
};


//- Read and add fields to the database
template<class FieldType>
void readAndAddFields(const fvMesh& mesh, const IOobjectList& objects)
{

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


//- Read and add all fields to the database
void readAndAddAllFields(const fvMesh& mesh)
{
    // Get all fields present at the current time
    IOobjectList objects(mesh, mesh.time().timeName());

    readAndAddFields<volScalarField>(mesh, objects);
    readAndAddFields<volVectorField>(mesh, objects);
    readAndAddFields<volSphericalTensorField>(mesh, objects);
    readAndAddFields<volSymmTensorField>(mesh, objects);
    readAndAddFields<volTensorField>(mesh, objects);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //- Add options
    timeSelector::addOptions(true, false);

    argList::addBoolOption
    (
        "updateAll",
        "Update size of all fields in the current time step"
    );
    argList::addBoolOption
    (
        "debug",
        "Output partial updates and additional fields"
    );
    argList::addBoolOption
    (
        "force3D",
        "Force 3-D refinement"
    );
    argList::addBoolOption
    (
        "noWrite",
        "Do not write fields"
    );
    argList::addBoolOption
    (
        "noZones",
        "Do not create zones"
    );
    argList::addBoolOption
    (
        "noRefine",
        "Do not refine"
    );

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Select time
    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const dictionary setFieldsDict(systemDict("setFieldsDict", args, mesh));

    // Update All fields in the current time folder
    // Resizes to make sure all fields are consistent with the mesh
    bool updateAll(args.optionFound("updateAll"));

    // Do not write fields
    // Usefull if a refined mesh is needed before mesh manipulation
    bool noWrite(args.optionFound("noWrite"));

    bool noZones(args.optionFound("noZones"));
    bool noRefine(args.optionFound("noRefine"));

    autoPtr<hexRef> meshCutter;
    bool refine = true;
    if (args.optionFound("force3D"))
    {
        meshCutter.set(new hexRef3D(mesh));
    }
    else if (mesh.nGeometricD() > 1)
    {
        meshCutter = hexRef::New(mesh);
    }
    else
    {
        // 1D so do not refine
        refine = false;
    }
    refine = refine && !noRefine;

    Switch debug(args.optionFound("debug"));

    label nBufferLayers(setFieldsDict.lookupOrDefault("nBufferLayers", 0));

    PackedBoolList protectedCell(mesh.nCells(), 0);
    if (refine)
    {
        initialize(mesh, protectedCell, meshCutter());
    }

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
    wordList fieldNames(setFieldsDict.lookupOrDefault("fields", wordList()));
    PtrList<volScalarField> fields(fieldNames.size());
    label fi = 0;
    forAll(fields, fieldi)
    {
        // Check the current time directory
        IOobject fieldHeader
        (
            fieldNames[fieldi],
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (fieldHeader.typeHeaderOk<volScalarField>(true))
        {
            fields.set
            (
                fi++,
                new volScalarField
                (
                    IOobject
                    (
                        fieldNames[fieldi],
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            WarningInFunction
                << "Field " << fieldNames[fieldi] << " specified for "
                << "setting refinement was not found. " << endl;
        }
    }
    fields.resize(fi);

    // Read in all fields to allow resizing
    if (updateAll)
    {
        readAndAddAllFields(mesh);
    }

    // Regions to refine based on a field
    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    //- List of sources (and backups if present)
    // Stored to reduce the number of reads
    PtrList<topoSetSource> sources(regions.size());
    PtrList<topoSetSource> backupSources(regions.size());
    HashPtrTable<topoSet> cellSets;
    HashPtrTable<topoSet> faceSets;
    HashPtrTable<topoSet> pointSets;

    labelList levels(regions.size(), 0);
    forAll(regions, regionI)
    {
        sources.set
        (
            regionI,
            topoSetSource::New
            (
                regions[regionI].keyword(),
                mesh,
                regions[regionI].dict()
            ).ptr()
        );
        if (regions[regionI].dict().found("backup"))
        {
            backupSources.set
            (
                regionI,
                topoSetSource::New
                (
                    regions[regionI].keyword(),
                    mesh,
                    regions[regionI].dict().subDict("backup")
                ).ptr()
            );
        }
        levels[regionI] =
            regions[regionI].dict().lookupOrDefault("level", 0);
    }

    label maxLevel =
    (
        levels.size() && !setFieldsDict.found("maxRefinement")
      ? max(levels)
      : setFieldsDict.lookupOrDefault<label>("maxRefinement", 0)
    );

    // Error fields is the same since it is looked up
    autoPtr<errorEstimator> EE;
    if (setFieldsDict.found("errorEstimator"))
    {
        EE = errorEstimator::New(mesh, setFieldsDict);
    }

    // Maximum number of iterations
    label iter = 0;
    label maxIter = 1;
    if (refine)
    {
        maxIter =
            max
            (
                1,
                max
                (
                    setFieldsDict.lookupOrDefault("maxIter", 2*maxLevel),
                    gMax(meshCutter->cellLevel())*2
                )
            );
    }

    // If debugging read in all to the necessary fields
    if (debug)
    {
        // Read default values and set fields
        if (setFieldsDict.found("defaultFieldValues"))
        {
            PtrList<topoSetList>
            (
                setFieldsDict.lookup("defaultFieldValues"),
                topoSetList::iNew(mesh, setFieldsDict)
            );
        }

        forAll(regions, regionI)
        {
            if (regions[regionI].dict().found("fieldValues"))
            {
                PtrList<topoSetList>
                (
                    regions[regionI].dict().lookup("fieldValues"),
                    topoSetList::iNew(mesh, setFieldsDict)
                );
            }
        }
    }

    // Flag for final iteration
    bool end = false;

    // Flag to initiate end
    bool prepareToStop = (maxIter == 1);

    while(!end)
    {
        if (debug)
        {
            runTime++;
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

        bool writeSets = !noZones && (end || debug);
        bool write = !noWrite && (end || debug);

        // Read default values and set fields
        if (setFieldsDict.found("defaultFieldValues"))
        {
            Info<< "Setting field default values" << endl;
            PtrList<topoSetList>
            (
                setFieldsDict.lookup("defaultFieldValues"),
                topoSetList::iNew
                (
                    mesh,
                    setFieldsDict,
                    identity(mesh.nCells()),
                    identity(mesh.nFaces()),
                    identity(mesh.nPoints()),
                    write
                )
            );

            Info<< endl;
        }

        // List of saved cells (per cell set)
        labelListList savedCells;

        Info<< "Setting field region values" << endl;
        forAll(regions, regionI)
        {
            const dictionary& regionDict =  regions[regionI].dict();
            labelList selectedCells;
            labelList selectedFaces;
            labelList selectedPoints;

            // Cell sets
            if (sources[regionI].setType() == topoSetSource::CELLSETSOURCE)
            {
                cellSet selectedCellSet
                (
                    mesh,
                    "cellSet",
                    mesh.nCells()/10+1  // Reasonable size estimate.
                );
                sources[regionI].applyToSet
                (
                    topoSetSource::NEW,
                    selectedCellSet
                );
                selectedCells = selectedCellSet.toc();

                if
                (
                    !returnReduce(selectedCells.size(), sumOp<label>())
                 && !end
                 && backupSources.set(regionI)
                )
                {
                    Info<< nl
                        << "    Expanding refinement region" << endl;
                    cellSet backupCellSet
                    (
                        mesh,
                        "cellSet",
                        mesh.nCells()/10+1  // Reasonable size estimate.
                    );

                    backupSources[regionI].applyToSet
                    (
                        topoSetSource::NEW,
                        backupCellSet
                    );
                    selectedCells = backupCellSet.toc();
                }

                if (!returnReduce(selectedCells.size(), sumOp<label>()))
                {
                    WarningInFunction
                        << "No cells were selected for using " << nl
                        << regions[regionI].name()
                        << regionDict
                        << "To expand searchable region add backup " << nl
                        << "Region or expand backup region." << endl;
                }
                selectedFaces = cellsToFaces(mesh, selectedCells);
                selectedPoints = cellsToPoints(mesh, selectedCells);
            }
            // Face sets
            else if
            (
                sources[regionI].setType() == topoSetSource::FACESETSOURCE
            )
            {
                faceSet selectedFaceSet
                (
                    mesh,
                    "faceSet",
                    mesh.nFaces()/10+1
                );
                sources[regionI].applyToSet
                (
                    topoSetSource::NEW,
                    selectedFaceSet
                );
                selectedFaces = selectedFaceSet.toc();

                if
                (
                    !returnReduce(selectedFaces.size(), sumOp<label>())
                 && !end
                 && backupSources.set(regionI)
                )
                {
                    Info<< nl
                        << "    Expanding refinement region" << endl;
                    faceSet backupFaceSet
                    (
                        mesh,
                        "faceSet",
                        mesh.nFaces()/10+1  // Reasonable size estimate.
                    );

                    backupSources[regionI].applyToSet
                    (
                        topoSetSource::NEW,
                        backupFaceSet
                    );
                    selectedFaces = backupFaceSet.toc();
                }

                if (!returnReduce(selectedFaces.size(), sumOp<label>()))
                {
                    WarningInFunction
                        << "No faces were selected for using " << nl
                        << regions[regionI].name()
                        << regionDict
                        << "To expand searchable region add backup " << nl
                        << "Region or expand backup region." << endl;
                }
                selectedCells = facesToCells(mesh, selectedFaces);
                selectedPoints = facesToPoints(mesh, selectedFaces);
            }

            // Point sets
            else if
            (
                sources[regionI].setType() == topoSetSource::POINTSETSOURCE
            )
            {
                pointSet selectedPointSet
                (
                    mesh,
                    "pointSet",
                    (mesh.points().size())/10+1
                );
                sources[regionI].applyToSet
                (
                    topoSetSource::NEW,
                    selectedPointSet
                );
                selectedPoints = selectedPointSet.toc();

                if
                (
                    !returnReduce(selectedPoints.size(), sumOp<label>())
                 && !end
                 && backupSources.set(regionI)
                )
                {
                    Info<< nl
                        << "    Expanding refinement region" << endl;
                    faceSet backupPointSet
                    (
                        mesh,
                        "faceSet",
                        mesh.nFaces()/10+1  // Reasonable size estimate.
                    );

                    backupSources[regionI].applyToSet
                    (
                        topoSetSource::NEW,
                        backupPointSet
                    );
                    selectedPoints = backupPointSet.toc();
                }

                if (!returnReduce(selectedPoints.size(), sumOp<label>()))
                {
                    WarningInFunction
                        << "No points were selected for using " << nl
                        << regions[regionI].name()
                        << regionDict
                        << "To expand searchable region add backup " << nl
                        << "Region or expand backup region." << endl;
                }
                selectedCells = pointsToCells(mesh, selectedPoints);
                selectedFaces = pointsToFaces(mesh, selectedPoints);
            }

            savedCells.append(selectedCells);
            bool set =
                selectedCells.size()
             || selectedFaces.size()
             || selectedPoints.size();
            reduce(set, orOp<bool>());

            if (regionDict.found("fieldValues") && set)
            {
                Info<< "    Selected "
                    << selectedCells.size() << " cells, "
                    << selectedFaces.size() << " faces, "
                    << selectedPoints.size() << " points" << endl;

                // Print the volume of the cells set
                scalar V = 0.0;
                forAll(selectedCells, celli)
                {
                    V += mesh.V()[selectedCells[celli]];
                }
                Info<< "    Set volume of cell: "
                    << returnReduce(V, sumOp<scalar>()) << " m^3"
                    << endl;

                PtrList<topoSetList>
                (
                    regionDict.lookup("fieldValues"),
                    topoSetList::iNew
                    (
                        mesh,
                        setFieldsDict,
                        selectedCells,
                        selectedFaces,
                        selectedPoints,
                        end || debug
                    )
                );

                if (writeSets)
                {
                    if (regionDict.found("setName"))
                    {
                        word setName = regionDict.lookup<word>("setName");
                        if (!cellSets.found(setName))
                        {
                            cellSets.insert
                            (
                                setName,
                                topoSet::New
                                (
                                    "cellSet",
                                    mesh,
                                    setName,
                                    IOobject::NO_READ
                                ).ptr()
                            );
                        }
                        if (!faceSets.found(setName))
                        {
                            faceSets.insert
                            (
                                setName,
                                topoSet::New
                                (
                                    "faceSet",
                                    mesh,
                                    setName,
                                    IOobject::NO_READ
                                ).ptr()
                            );
                        }
                        if (!pointSets.found(setName))
                        {
                            pointSets.insert
                            (
                                setName,
                                topoSet::New
                                (
                                    "pointSet",
                                    mesh,
                                    setName,
                                    IOobject::NO_READ
                                ).ptr()
                            );
                        }
                        topoSet& cSet = *cellSets[setName];
                        topoSet& fSet = *faceSets[setName];
                        topoSet& pSet = *pointSets[setName];
                        forAll(selectedCells, i)
                        {
                            cSet.insert(selectedCells[i]);
                        }
                        forAll(selectedFaces, i)
                        {
                            fSet.insert(selectedFaces[i]);
                        }
                        forAll(selectedPoints, i)
                        {
                            pSet.insert(selectedPoints[i]);
                        }
                    }
                }
            }
            Info<< endl;
        }

        // Update boundary conditions of fields used for refinement
        forAll(fields, fieldi)
        {
            fields[fieldi].correctBoundaryConditions();
        }

        // Update error and mesh if not the final iteration
        if (!end && refine)
        {
            // Use error estimator to calculate error
            if (EE.valid())
            {
                EE->update();
            }
            // Check for any difference
            else
            {
                calcFaceDiff(error, fields);
            }
            labelList maxCellLevel(mesh.nCells(), -1);
            forAll(regions, regionI)
            {
                // Set keyword to lookup for refinement
                word refineKeyword;
                if (sources[regionI].setType() == topoSetSource::CELLSETSOURCE)
                {
                    refineKeyword = "refineInternal";
                }
                else if
                (
                    sources[regionI].setType() == topoSetSource::FACESETSOURCE
                )
                {
                    refineKeyword = "refineFaces";
                }
                else if
                (
                    sources[regionI].setType() == topoSetSource::POINTSETSOURCE
                )
                {
                    refineKeyword = "refinePoints";
                }

                // Do not use the max level, use current
                // Order is important in the definitions of regions
                Switch overwrite =
                    regions[regionI].dict().lookupOrDefault<Switch>
                    (
                        "overwriteLevel",
                        false
                    );

                // Set actual max cell level
                forAll(savedCells[regionI], celli)
                {
                    if (overwrite)
                    {
                        maxCellLevel[savedCells[regionI][celli]] = levels[regionI];
                    }
                    else
                    {
                        maxCellLevel[savedCells[regionI][celli]] =
                            max
                            (
                                maxCellLevel[savedCells[regionI][celli]],
                                levels[regionI]
                            );
                    }
                }

                // Set specified cells to be refined
                if
                (
                    regions[regionI].dict().lookupOrDefaultBackwardsCompatible<Switch>
                    (
                        {refineKeyword, "refine"},
                        false
                    ) && !EE.valid()
                )
                {
                    forAll(savedCells[regionI], celli)
                    {
                        error[savedCells[regionI][celli]] = 1.0;
                    }
                }

                // Extend refinement by nBufferLayers
                for (label i = 0; i < nBufferLayers; i++)
                {
                    extendMaxCellLevel
                    (
                        mesh,
                        savedCells[regionI],
                        maxCellLevel,
                        levels[regionI],
                        i == 0
                    );
                }
            }

            labelList maxRefinement(mesh.nCells(), maxLevel);
            if (EE.valid())
            {
                maxRefinement = EE->maxRefinement();
            }

            // Set the maxCell level
            forAll(maxCellLevel, celli)
            {
                if (maxCellLevel[celli] < 0)
                {
                    maxCellLevel[celli] = maxRefinement[celli];
                }
            }

            // Mark cells greater than the max cell level for unrefinment
            forAll(error, celli)
            {
                if (meshCutter->cellLevel()[celli] > maxCellLevel[celli])
                {
                    error[celli] = -1.0;
                }
            }

            // Write fields and mesh if using debug
            if (debug)
            {
                bool writeOk = (mesh.write() && meshCutter->write());
                volScalarField scalarCellLevel
                (
                    IOobject
                    (
                        "cellLevel",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false
                    ),
                    mesh,
                    dimensionedScalar(dimless, 0)
                );
                volScalarField scalarMaxCellLevel
                (
                    IOobject
                    (
                        "maxCellLevel",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false
                    ),
                    mesh,
                    dimensionedScalar(dimless, 0)
                );

                const labelList& cellLevel = meshCutter->cellLevel();

                forAll(cellLevel, celli)
                {
                    scalarCellLevel[celli] = cellLevel[celli];
                    scalarMaxCellLevel[celli] = maxCellLevel[celli];
                }
                writeOk = writeOk && scalarCellLevel.write();
                error.write();
                scalarMaxCellLevel.write();
            }

            // Update mesh (return if mesh changes)
            Info<< nl << "Updating mesh" << endl;
            prepareToStop =
               !update
                (
                    mesh,
                    error,
                    setFieldsDict,
                    meshCutter(),
                    protectedCell,
                    maxCellLevel
                );
            // Force refinement data to go to the current time directory.
            meshCutter->setInstance(runTime.timeName());
        }
        iter++;
    }

    // Write all fields
    if (updateAll)
    {
        runTime.write();
    }

    if ((cellSets.size() || faceSets.size()) && !noZones)
    {
        List<cellZone*> cellZones(cellSets.size());
        List<faceZone*> faceZones(faceSets.size());
        List<pointZone*> pointZones(pointSets.size());
        label i = 0;
        forAllConstIter
        (
            HashPtrTable<topoSet>,
            cellSets,
            iter
        )
        {
            cellZones[i] = new cellZone
            (
                iter()->name(),
                iter()->toc(),
                i,
                mesh.cellZones()
            );
        }

        i = 0;
        forAllConstIter
        (
            HashPtrTable<topoSet>,
            faceSets,
            iter
        )
        {
            faceZones[i] = new faceZone
            (
                iter()->name(),
                iter()->toc(),
                boolList(iter()->size(), false),
                i,
                mesh.faceZones()
            );
            i++;
        }

        i = 0;
        forAllConstIter
        (
            HashPtrTable<topoSet>,
            pointSets,
            iter
        )
        {
            pointZones[i] = new pointZone
            (
                iter()->name(),
                iter()->toc(),
                i,
                mesh.pointZones()
            );
            i++;
        }
        mesh.pointZones().clear();
        mesh.faceZones().clear();
        mesh.cellZones().clear();

        mesh.addZones(pointZones, faceZones, cellZones);
    }

    if (refine)
    {
        // Handle cell level (as volScalarField) explicitly
        // Important when using mapFields or rotateFields
        const labelIOList& cellLevel = meshCutter->cellLevel();
        if (mesh.foundObject<volScalarField>("cellLevel"))
        {
            volScalarField& scalarCellLevel =
                mesh.lookupObjectRef<volScalarField>("cellLevel");
            forAll(cellLevel, celli)
            {
                scalarCellLevel[celli] = cellLevel[celli];
            }
            scalarCellLevel.write();
        }
        else if (debug)
        {
            volScalarField scalarCellLevel
            (
                IOobject
                (
                    "cellLevel",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimless, 0)
            );
            forAll(cellLevel, celli)
            {
                scalarCellLevel[celli] = cellLevel[celli];
            }
            scalarCellLevel.write();
        }

        // Write mesh and cell levels
        meshCutter->write();
        mesh.write();

        //- Write points0 field to time directory
        pointIOField points0
        (
            IOobject
            (
                "points0",
                runTime.timeName(),
                polyMesh::meshSubDir,
                mesh
            ),
            mesh.points()
        );
        points0.write();

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
