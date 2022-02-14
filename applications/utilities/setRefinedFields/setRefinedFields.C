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
#include "backupTopoSetSource.H"
#include "topoSetList.H"
#include "volFields.H"
#include "systemDict.H"

#include "fvMeshHexRefiner.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"
#include "errorEstimator.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "calcAngleFraction.H"
#include "IOobjectList.H"
#include "FieldSetType.H"

using namespace Foam;


void calcFaceDiff
(
    volScalarField& error,
    const PtrList<volScalarField>& fields
)
{
    volScalarField errorOrig(error);
    const labelUList& owner = error.mesh().owner();
    const labelUList& neighbour = error.mesh().neighbour();
    const label nInternalFaces = error.mesh().nInternalFaces();

    forAll(fields, fieldi)
    {
        const volScalarField& field = fields[fieldi];
        for (label facei = 0; facei < nInternalFaces; facei++)
        {
            label own = owner[facei];
            label nei = neighbour[facei];
            scalar e =
                pos0(mag(field[own] - field[nei]) - 1e-6)
              - neg(mag(field[own] - field[nei]) - 1e-6);

            error[own] = max(error[own], e);
            error[nei] = max(error[nei], e);
        }

        // Boundary faces
        forAll(error.boundaryField(), patchi)
        {
            if (error.boundaryField()[patchi].coupled())
            {
                const fvPatch& p = field.boundaryField()[patchi].patch();

                const labelUList& faceCells = p.faceCells();
                scalarField fp
                (
                    field.boundaryField()[patchi].patchInternalField()
                );

                scalarField fn
                (
                    field.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(faceCells, facei)
                {
                    scalar e =
                        pos0(mag(fp[facei] - fn[facei]) - 1e-6)
                      - neg(mag(fp[facei] - fn[facei]) - 1e-6);
                    error[faceCells[facei]] =
                        max(error[faceCells[facei]], e);
                }
            }
        }
    }
}


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
    FatalErrorInFunction
        << "Could not determine geometry type from " << type << endl
        << abort(FatalError);
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
    FatalErrorInFunction
        << "Could not determine primitive type from " << type << endl
        << abort(FatalError);
    return UNKNOWN_PRIM;
}


class fieldSetList
{

public:

    fieldSetList()
    {}

    autoPtr<fieldSetList> clone() const
    {
        return autoPtr<fieldSetList>(new fieldSetList());
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
            const word fieldName(is);
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
            if (!returnReduce(elms.size(), sumOp<label>()) && !force_)
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

        autoPtr<fieldSetList> operator()(Istream& is) const
        {
            word fieldSetType(is);
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

            return autoPtr<fieldSetList>(new fieldSetList());
        }
    };
};


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
        "forceHex8",
        "Force use of standard OpenFOAM hexRef8 refinement"
    );
    argList::addBoolOption
    (
        "noWrite",
        "Do not write fields"
    );
    argList::addBoolOption
    (
        "noFields",
        "Do set fields"
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

    const dictionary setFieldsDict(systemDict("setFieldsDict", args, mesh));

    // Update All fields in the current time folder
    // Resizes to make sure all fields are consistent with the mesh
    bool updateAll(args.optionFound("updateAll"));

    // Do not write fields
    // Usefull if a refined mesh is needed before mesh manipulation
    bool noFields(args.optionFound("noFields"));

    // Do not write fields
    // Usefull if a refined mesh is needed before mesh manipulation
    bool noWrite(args.optionFound("noWrite") || noFields);

    bool noRefine(args.optionFound("noRefine"));
    bool overwrite(args.optionFound("overwrite"));
    bool noHistory(args.optionFound("noHistory"));
    bool debug(args.optionFound("debug"));

    if (overwrite && debug)
    {
        FatalErrorInFunction
            << "Cannot overwrite mesh in debug mode" << endl
            << abort(FatalError);
    }
    if (overwrite)
    {
        mesh.setInstance(runTime.constant());
    }

    //- Is the mesh balanced
    bool balance = false;
    autoPtr<fvMeshRefiner> refiner;
    if (!noRefine)
    {
        dictionary refineDict
        (
            setFieldsDict.optionalSubDict
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
        if (Pstream::parRun())
        {
            balance = refineDict.lookupOrDefault("balance", true);
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

    wordList fieldNames;
    if (!noFields)
    {
        fieldNames = setFieldsDict.lookupOrDefault("fields", wordList());
    }
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldHeader.typeHeaderOk<volScalarField>(true))
        {
            fields.set
            (
                fi++,
                new volScalarField(fieldHeader, mesh)
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
    const scalar angleFraction = calcAngleFraction(mesh);

    // Read in all fields to allow resizing
    if (updateAll)
    {
        readAndAddAllFields(mesh);
    }

    //- List of sources (and backups if present)
    // Stored to reduce the number of reads
    PtrList<backupTopoSetSource> regions
    (
        setFieldsDict.lookup("regions"),
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
                setFieldsDict.lookupOrDefault
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

    // If debugging read in all to the necessary fields
    if (debug && !noFields)
    {
        // Read default values and set fields
        if (setFieldsDict.found("defaultFieldValues"))
        {
            PtrList<fieldSetList>
            (
                setFieldsDict.lookup("defaultFieldValues"),
                fieldSetList::iNew(mesh, setFieldsDict)
            );
        }

        forAll(regions, regionI)
        {
            if (regions[regionI].dict().found("fieldValues"))
            {
                PtrList<fieldSetList>
                (
                    regions[regionI].dict().lookup("fieldValues"),
                    fieldSetList::iNew(mesh, setFieldsDict)
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
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        if (maxIter <= iter)
        {
            prepareToStop = true;
            WarningInFunction << "Reach maximum number of iteration" << endl;
        }

        error = -1.0;

        // Check if this is the final iteration so correct cell shapes are used
        if (prepareToStop)
        {
            end = true;
        }

        bool write = !noWrite && (end || debug);

        // Read default values and set fields
        if (setFieldsDict.found("defaultFieldValues") && !noFields)
        {
            Info<< "Setting field default values" << endl;
            PtrList<fieldSetList>
            (
                setFieldsDict.lookup("defaultFieldValues"),
                fieldSetList::iNew
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
        labelListList savedCells(regions.size());
        labelListList savedFaces(regions.size());
        boolListList savedFlipMaps(regions.size());
        labelListList savedPoints(regions.size());

        Info<< "Setting field region values" << endl;
        forAll(regions, regionI)
        {
            const dictionary& regionDict =  regions[regionI].dict();
            regions[regionI].allowBackup(!end);
            regions[regionI].updateSets();

            const labelList& selectedCells = regions[regionI].selectedCells();
            const labelList& selectedFaces = regions[regionI].selectedFaces();
            const boolList& selectedFlipMaps = regions[regionI].flipMap();
            const labelList& selectedPoints = regions[regionI].selectedPoints();

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

            savedCells[regionI] = selectedCells;
            savedFaces[regionI] = selectedFaces;
            savedFlipMaps[regionI] = selectedFlipMaps;
            savedPoints[regionI] = selectedPoints;

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

                // Print the volume of the cells set
                scalar V = 0.0;
                forAll(selectedCells, celli)
                {
                    V += mesh.V()[selectedCells[celli]];
                }
                reduce(V, sumOp<scalar>());
                Info<< "    Set volume of cell: "
                    << V << " m^3";
                if (angleFraction != 1.0)
                {
                    Info<< " (" << V/angleFraction << " actual)";
                }
                Info<< endl;

                if (regionDict.found("fieldValues") && !noFields)
                {
                    PtrList<fieldSetList>
                    (
                        regionDict.lookup("fieldValues"),
                        fieldSetList::iNew
                        (
                            mesh,
                            setFieldsDict,
                            selectedCells,
                            selectedFaces,
                            selectedPoints,
                            end || debug
                        )
                    );
                }
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

        // Update boundary conditions of fields used for refinement
        forAll(fields, fieldi)
        {
            fields[fieldi].correctBoundaryConditions();
        }

        // Update error and mesh if not the final iteration
        if (refine)
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

                if (!EE.valid())
                {
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
                 && scalarMaxCellLevel.write()
                 && error.write();

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
    topoSets.transferZones(true);

    bool writeMesh = topoSets.writeSets();

    if (refine && !debug)
    {
        // Write mesh and cell levels
        if (overwrite)
        {
            mesh.setInstance(runTime.constant());
        }

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

    if (noHistory)
    {
        refiner.clear();
    }

    // Write all fields
    if (updateAll)
    {
        runTime.write();
    }
    if (writeMesh)
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
