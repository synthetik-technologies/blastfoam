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
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "volFields.H"

#include "hexRef.H"
#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"

using namespace Foam;

#include "refineFunctions.H"

template<class Type>
bool setCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream,
    const bool write
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    bool useOrig = false;
    if (fieldTypeDesc == fieldType::typeName + "InitialValue")
    {
        useOrig = true;
    }
    else if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);
    word fieldOrigName(fieldName + "Orig");
    fieldType* field = 0;
    if (mesh.foundObject<fieldType>(fieldName))
    {
        field = &mesh.lookupObjectRef<fieldType>(fieldName);
    }
    else
    {
        // Check the current time directory
        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check the "constant" directory
        if (!fieldHeader.typeHeaderOk<fieldType>(true))
        {
            fieldHeader = IOobject
            (
                fieldName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            );
        }

        // Check field exists
        if (fieldHeader.typeHeaderOk<fieldType>(true))
        {
            field = new fieldType(fieldHeader, mesh);
        }
    }
    if (field && useOrig)
    {
        Info<< "    Setting internal values of "
            << field->headerClassName()
            << " " << fieldName << " to initial values" << endl;
    }
    else if (field)
    {
        Info<< "    Setting internal values of "
            << field->headerClassName()
            << " " << fieldName << endl;
    }
    else
    {
        WarningInFunction
            << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    if (useOrig)
    {
        fieldType* origFieldPtr = 0;
        if (mesh.foundObject<fieldType>(fieldOrigName))
        {
            origFieldPtr = &mesh.lookupObjectRef<fieldType>(fieldOrigName);
        }
        else
        {
            origFieldPtr = new fieldType(fieldOrigName, *field);
            origFieldPtr->store(origFieldPtr);
        }
        (*field) = (*origFieldPtr);
        return true;
    }

    const Type& value = pTraits<Type>(fieldValueStream);

    if (selectedCells.size() == field->size())
    {
        field->primitiveFieldRef() = value;
    }
    else
    {
        forAll(selectedCells, celli)
        {
            (*field)[selectedCells[celli]] = value;
        }
    }

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& fieldBf = field->boundaryFieldRef();

    forAll(field->boundaryField(), patchi)
    {
        fieldBf[patchi] = fieldBf[patchi].patchInternalField();
    }

    if (field->name() != "error" && write)
    {
        if (!field->write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field->name() << exit(FatalError);
        }
    }

    return true;
}


class setCellField
{

public:

    setCellField()
    {}

    autoPtr<setCellField> clone() const
    {
        return autoPtr<setCellField>(new setCellField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;
        const bool write_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells, const bool write)
        :
            mesh_(mesh),
            selectedCells_(selectedCells),
            write_(write)
        {}

        autoPtr<setCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setCellFieldType<scalar>
                        (fieldType, mesh_, selectedCells_, fieldValues, write_)
                 || setCellFieldType<vector>
                        (fieldType, mesh_, selectedCells_, fieldValues, write_)
                 || setCellFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues, write_)
                 || setCellFieldType<symmTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues, write_)
                 || setCellFieldType<tensor>
                        (fieldType, mesh_, selectedCells_, fieldValues, write_)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setCellField>(new setCellField());
        }
    };
};


template<class Type>
bool setFaceFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedFaces,
    Istream& fieldValueStream,
    const bool write
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);
    fieldType* field = 0;
    if (mesh.foundObject<fieldType>(fieldName))
    {
        field =  &mesh.lookupObjectRef<fieldType>(fieldName);
    }
    else
    {
        // Check the current time directory
        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check the "constant" directory
        if (!fieldHeader.typeHeaderOk<fieldType>(true))
        {
            fieldHeader = IOobject
            (
                fieldName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            );
        }
    }
    if (field)
    {
        Info<< "    Setting internal values of "
            << field->headerClassName()
            << " " << fieldName << endl;
    }
    else
    {
        WarningInFunction
            << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    Info<< "    Setting patchField values of "
        << field->headerClassName()
        << " " << fieldName << endl;

    const Type& value = pTraits<Type>(fieldValueStream);

    // Create flat list of selected faces and their value.
    Field<Type> allBoundaryValues(mesh.nFaces()-mesh.nInternalFaces());
    forAll(field->boundaryField(), patchi)
    {
        SubField<Type>
        (
            allBoundaryValues,
            field->boundaryField()[patchi].size(),
            field->boundaryField()[patchi].patch().start()
          - mesh.nInternalFaces()
        ) = field->boundaryField()[patchi];
    }

    // Override
    bool hasWarned = false;
    labelList nChanged
    (
        returnReduce(field->boundaryField().size(), maxOp<label>()),
        0
    );
    forAll(selectedFaces, i)
    {
        label facei = selectedFaces[i];
        if (mesh.isInternalFace(facei))
        {
            if (!hasWarned)
            {
                hasWarned = true;
                WarningInFunction
                    << "Ignoring internal face " << facei
                    << ". Suppressing further warnings." << endl;
            }
        }
        else
        {
            label bFacei = facei-mesh.nInternalFaces();
            allBoundaryValues[bFacei] = value;
            nChanged[mesh.boundaryMesh().patchID()[bFacei]]++;
        }
    }

    Pstream::listCombineGather(nChanged, plusEqOp<label>());
    Pstream::listCombineScatter(nChanged);

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& fieldBf = field->boundaryFieldRef();

    // Reassign.
    forAll(field->boundaryField(), patchi)
    {
        if (nChanged[patchi] > 0)
        {
            Info<< "    On patch "
                << field->boundaryField()[patchi].patch().name()
                << " set " << nChanged[patchi] << " values" << endl;
            fieldBf[patchi] == SubField<Type>
            (
                allBoundaryValues,
                fieldBf[patchi].size(),
                fieldBf[patchi].patch().start()
              - mesh.nInternalFaces()
            );
        }
    }

    if (field->name() != "error" && write)
    {
        if (!field->write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field->name() << exit(FatalError);
        }
    }

    return true;
}


class setFaceField
{

public:

    setFaceField()
    {}

    autoPtr<setFaceField> clone() const
    {
        return autoPtr<setFaceField>(new setFaceField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedFaces_;
        const bool write_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedFaces, const bool write)
        :
            mesh_(mesh),
            selectedFaces_(selectedFaces),
            write_(write)
        {}

        autoPtr<setFaceField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setFaceFieldType<scalar>
                        (fieldType, mesh_, selectedFaces_, fieldValues, write_)
                 || setFaceFieldType<vector>
                        (fieldType, mesh_, selectedFaces_, fieldValues, write_)
                 || setFaceFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues, write_)
                 || setFaceFieldType<symmTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues, write_)
                 || setFaceFieldType<tensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues, write_)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setFaceField>(new setFaceField());
        }
    };
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("setFieldsDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary setFieldsDict(dictIO);

    autoPtr<hexRef> meshCutter;
    if (setFieldsDict.lookupOrDefault("force3D", Switch(false)))
    {
        meshCutter.set(new hexRef8(mesh));
    }
    else
    {
        meshCutter = hexRef::New(mesh);
    }
    PackedBoolList protectedCell(mesh.nCells(), 0);
    initialize(mesh, protectedCell, meshCutter());
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
    wordList fieldNames(setFieldsDict.lookup("fields"));
    PtrList<volScalarField> fields(fieldNames.size());
    forAll(fields, fieldi)
    {
        fields.set
        (
            fieldi,
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

    // Regions to refine based on a field
    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    //- List of source
    PtrList<topoSetSource> sources(regions.size());
    PtrList<topoSetSource> backupSources(regions.size());
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
            regions[regionI].dict().lookupType<label>("level");
    }
    label maxRefinement(max(levels));

    bool end = false;
    bool prepareToStop = false;
    labelList nOldCells(regions.size(), -1);

    label iter = 0;
    label maxIter = setFieldsDict.lookupOrDefault("maxIter", 2*maxRefinement);
    while(!end)
    {
        if (maxIter <= iter)
        {
            prepareToStop = true;
        }

        error = -10.0;
        if (prepareToStop)
        {
            end = true;
        }
        if (setFieldsDict.found("defaultFieldValues"))
        {
            Info<< "Setting field default values" << endl;
            PtrList<setCellField> defaultFieldValues
            (
                setFieldsDict.lookup("defaultFieldValues"),
                setCellField::iNew(mesh, labelList(mesh.nCells()), end)
            );
            Info<< endl;
        }

        labelListList savedCells;

        Info<< "Setting field region values" << endl;
        forAll(regions, regionI)
        {
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

                labelList cells = selectedCellSet.toc();

                if
                (
                    returnReduce
                    (
                        selectedCellSet.toc().size(), sumOp<label>()
                    ) > 0
                )
                {
                    cells = selectedCellSet.toc();
                }
                else if (!end && backupSources.set(regionI))
                {
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
                    cells = backupCellSet.toc();
                }

                //- Mark for possible refinement
                {
                    const labelUList& owner = error.mesh().owner();
                    const labelUList& neighbour = error.mesh().neighbour();

                    forAll(cells, celli)
                    {
                        const labelList& faces
                        (
                            error.mesh().cells()[cells[celli]]
                        );
                        forAll(faces, facei)
                        {
                            if (faces[facei] < error.mesh().nInternalFaces())
                            {
                                error[owner[faces[facei]]] = 0.1;
                                error[neighbour[faces[facei]]] = 0.1;
                            }
                        }
                    }
                }

                if (regions[regionI].dict().found("fieldValues"))
                {
                    PtrList<setCellField> fieldValues
                    (
                        regions[regionI].dict().lookup("fieldValues"),
                        setCellField::iNew(mesh, cells, end)
                    );
                }
                nOldCells[regionI] = cells.size();
                savedCells.append(cells);

            }
            else if (sources[regionI].setType() == topoSetSource::FACESETSOURCE)
            {
                faceSet selectedFaceSet
                (
                    mesh,
                    "faceSet",
                    (mesh.nFaces()-mesh.nInternalFaces())/10+1
                );

                sources[regionI].applyToSet
                (
                    topoSetSource::NEW,
                    selectedFaceSet
                );

                labelList selectedFaces(selectedFaceSet.toc());
                if (regions[regionI].dict().found("fieldValues"))
                {
                    PtrList<setFaceField> fieldValues
                    (
                        regions[regionI].dict().lookup("fieldValues"),
                        setFaceField::iNew(mesh, selectedFaces, end)
                    );
                }

                labelList faceCells(selectedFaces.size()*2);
                label nFaces = 0;
                forAll(selectedFaces, i)
                {
                    label facei = selectedFaces[i];
                    faceCells[nFaces] = mesh.faceOwner()[facei];
                    nFaces++;
                    if (mesh.isInternalFace(facei))
                    {
                        faceCells[nFaces] = mesh.faceNeighbour()[facei];
                        nFaces++;
                    }
                }
                faceCells.resize(nFaces);
                savedCells.append(faceCells);
            }
        }
        forAll(fields, fieldi)
        {
            fields[fieldi].correctBoundaryConditions();
        }

        if (!end)
        {
            // Refine internal cells
            calcFaceDiff(error, fields);
            labelList maxCellLevel(mesh.nCells(), -1);
            forAll(regions, regionI)
            {
                word refineKeyword;
                if (sources[regionI].setType() == topoSetSource::CELLSETSOURCE)
                {
                    refineKeyword = "refineInternal";
                }
                else
                {
                    refineKeyword = "refineFaces";
                }
                forAll(savedCells[regionI], celli)
                {
                    maxCellLevel[savedCells[regionI][celli]] =
                        max
                        (
                            maxCellLevel[savedCells[regionI][celli]],
                            levels[regionI]
                        );
                }

                if
                (
                    regions[regionI].dict().lookupOrDefault<Switch>
                    (
                        refineKeyword,
                        false
                    )
                )
                {
                    forAll(savedCells[regionI], celli)
                    {
                        error[savedCells[regionI][celli]] = 1.0;
                    }
                }


                label nBufferLayers
                (
                    setFieldsDict.lookupType<label>("nBufferLayers")
                );
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
            forAll(maxCellLevel, celli)
            {
                if (maxCellLevel[celli] < 0)
                {
                    maxCellLevel[celli] = maxRefinement;
                }
            }

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

    const labelList& cellLevel = meshCutter->cellLevel();

    forAll(cellLevel, celli)
    {
        scalarCellLevel[celli] = cellLevel[celli];
    }
    writeOk = writeOk && scalarCellLevel.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
