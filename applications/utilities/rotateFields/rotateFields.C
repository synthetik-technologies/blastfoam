/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
10-29-2020 Jeff Heylmun     | Specialized mapFields to rotate axisymmetric cases
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

Application
    rotateFields

Description
    Rotate mesh and fields from 1-D to 2-D or 2-D to 3-D. Only for
    axisymmetric cases. Mapping from parallel cases is not currently supported

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "labelVector.H"
#include "wedgeFvPatch.H"
#include "IOobjectList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
wordList createBoundaryTypes
(
    const GeometricField<Type, fvPatchField, volMesh>& src,
    const fvMesh& targetMesh
)
{
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        srcBoundary = src.boundaryField();

    HashTable<word> srcBoundaryTypes;
    forAll(srcBoundary, patchi)
    {
        srcBoundaryTypes.insert
        (
            srcBoundary[patchi].patch().name(),
            srcBoundary[patchi].type()
        );
    }

    wordList targetBoundaryTypes
    (
        targetMesh.boundary().size(),
        calculatedFvPatchField<Type>::typeName
    );
    forAll(targetMesh.boundary(), patchi)
    {
        if (srcBoundaryTypes.found(targetMesh.boundaryMesh()[patchi].name()))
        {
            targetBoundaryTypes[patchi] =
                srcBoundaryTypes[targetMesh.boundaryMesh()[patchi].name()];
        }
    }
    return targetBoundaryTypes;
}


template<class Type>
void mapVolFields
(
    const fvMesh& sourceMesh,
    const fvMesh& targetMesh,
    const labelList& cellMap,
    const IOobjectList& objects,
    const tensorField& R,
    const HashSet<word>& mapFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    IOobjectList fields = objects.lookupClass(fieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            targetMesh.time().timeName(),
            targetMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        bool map = false;
        fieldType fieldSource(*fieldIter(), sourceMesh);
        autoPtr<fieldType> fieldTargetPtr;
        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            map = true;
            fieldTargetPtr.set
            (
                new fieldType
                (
                    fieldTargetIOobject,
                    targetMesh
                )
            );
        }
        else if (mapFields.found(fieldTargetIOobject.name()))
        {
            map = true;
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;
            fieldTargetPtr.set
            (
                new fieldType
                (
                    fieldTargetIOobject,
                    targetMesh,
                    dimensioned<Type>
                    (
                        "0",
                        fieldSource.dimensions(),
                        pTraits<Type>::zero
                    ),
                    createBoundaryTypes<Type>(fieldSource, targetMesh)
                )
            );
        }

        if (map)
        {
            Info<< "    mapping " << fieldIter()->name() << endl;

            // Read fieldTarget
            fieldType& fieldTarget = fieldTargetPtr();

            forAll(cellMap, celli)
            {
                label cellj = cellMap[celli];
                if (cellj != -1)
                {
                    Type v = fieldSource[cellj];
                    fieldTarget[celli] = (R[celli] & v);
                }
            }
            fieldTarget.write();
        }
    }
}

void mapVolScalarFields
(
    const fvMesh& sourceMesh,
    const fvMesh& targetMesh,
    const labelList& cellMap,
    const IOobjectList& objects,
    const HashSet<word>& mapFields
)
{
    typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
    IOobjectList fields = objects.lookupClass(fieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            targetMesh.time().timeName(),
            targetMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // Read field fieldSource
        fieldType fieldSource(*fieldIter(), sourceMesh);
        autoPtr<fieldType> fieldTargetPtr;

        bool map = false;
        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            map = true;
            fieldTargetPtr.set
            (
                new fieldType
                (
                    fieldTargetIOobject,
                    targetMesh
                )
            );
        }
        else if (mapFields.found(fieldTargetIOobject.name()))
        {
            map = true;
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;
            fieldTargetPtr.set
            (
                new fieldType
                (
                    fieldTargetIOobject,
                    targetMesh,
                    dimensioned<scalar>
                    (
                        "0",
                        fieldSource.dimensions(),
                        0.0
                    ),
                    createBoundaryTypes<scalar>(fieldSource, targetMesh)
                )
            );
        }
        if (map)
        {
            Info<< "    mapping " << fieldIter()->name() << endl;

            // Read fieldTarget
            fieldType& fieldTarget = fieldTargetPtr();
            forAll(cellMap, celli)
            {
                label cellj = cellMap[celli];
                if (cellj != -1)
                {
                    fieldTarget[celli] = fieldSource[cellj];
                }
            }
            fieldTarget.correctBoundaryConditions();
            fieldTarget.write();
        }
    }
}

void mapFields
(
    const fvMesh& sourceMesh,
    const fvMesh& targetMesh,
    const labelList& cellMap,
    const tensorField& R,
    const HashSet<word>& additionalFields
)
{
    IOobjectList objects(sourceMesh, sourceMesh.time().timeName());

    mapVolScalarFields(sourceMesh, targetMesh, cellMap, objects, additionalFields);
    mapVolFields<vector>(sourceMesh, targetMesh, cellMap, objects, R, additionalFields);
//     mapVolFields<sphericalTensor>(sourceMesh, targetMesh, cellMap, objects, R);
//     mapVolFields<symmTensor>(sourceMesh, targetMesh, cellMap, objects, R);
    mapVolFields<tensor>(sourceMesh, targetMesh, cellMap, objects, R, additionalFields);
}


Foam::Pair<Foam::vector> calculateAxis(const fvMesh& mesh)
{
    Pair<vector> axis(vector::one, Zero);
    vector& rAxis = axis[0];
    vector& yAxis = axis[1];

    List<vector> foundAxis;
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            const wedgePolyPatch& wp = dynamicCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchi]
            );
            bool found = false;
            vector a = cmptMag(wp.axis());
            forAll(foundAxis, ai)
            {
                if (mag(a - foundAxis[ai]) < 1e-10)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                foundAxis.append(a);
                rAxis -= cmptMag(wp.centreNormal()) + a;
                yAxis += a;
            }
        }
    }
    forAll(yAxis, cmpti)
    {
        rAxis[cmpti] = min(rAxis[cmpti]*pos(rAxis[cmpti] - small), 1.0);
        yAxis[cmpti] = min(yAxis[cmpti]*pos(yAxis[cmpti] - small), 1.0);
    }
    return axis;
}


Foam::vector calculateCentre(const fvMesh& mesh)
{
    return
        gSum
        (
            cmptMultiply
            (
                mesh.C().primitiveField()*mesh.V().field(),
                calculateAxis(mesh)[1]
            )
        )/gSum(mesh.V().field());
}


void calcMapAndR
(
    const fvMesh& sourceMesh,
    const fvMesh& targetMesh,
    const scalar& maxR,
    const vector& sourceCentre,
    const vector& targetCentre,
    const vector& rotationAxis,
    const vector& rAxis,
    labelList& cellMap,
    tensorField& R
)
{
    label nSourceD = sourceMesh.nGeometricD();
    vector sourceD(sourceMesh.geometricD());
    vector targetD(targetMesh.geometricD());

    // Create map
    cellMap = labelList(targetMesh.nCells(), -1);
    R = tensorField(targetMesh.nCells(), tensor::I);
    forAll(cellMap, celli)
    {
        // Get the position on the target mesh
        vector ptTarget =
            cmptMultiply
            (
                targetMesh.cellCentres()[celli],
                targetD
            );

        // Offset from the center of the source mesh
        vector nTarget = (ptTarget - targetCentre);

        // Radius
        scalar r = mag(nTarget);

        // Remove points outside of the maximum radius
        if (r > maxR)
        {
            r = great;
        }

        // Offset from the source mesh center (only solved directions)
        vector nSource(Zero);
        if (nSourceD == 1)
        {
            nSource = r*rAxis;
        }
        else
        {
            scalar y = nTarget & rotationAxis;
            scalar x = mag((nTarget - y*rotationAxis));
            nSource = y*rotationAxis + x*rAxis;
        }

        // Actual point on the source mesh
        vector ptSource = nSource + sourceCentre;

        // Map from the source mesh to the target mesh
        cellMap[celli] = sourceMesh.findCell(ptSource);

        // Extend radius is the target point is outside of the source mesh
        if (maxR < 0 && cellMap[celli] < 0)
        {
            cellMap[celli] = sourceMesh.findNearestCell(ptSource);
        }

        // Calculate the rotation matrix for vectors and tensors
        if (cellMap[celli] >= 0)
        {
            vector a = nSource/mag(nSource);
            vector b = nTarget/mag(nTarget);
            vector v = (a ^ b);
            scalar c = a & b;
            tensor A
            (
                0.0, -v[2], v[1],
                v[2], 0.0, -v[0],
                -v[1], v[0], 0.0
            );
            R[celli] = tensor::I + A + (A & A)/(1.0 + c);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //- Add options
    timeSelector::addOptions(true, false);

    argList::addNote
    (
        "map volume fields from one mesh to another"
    );
    argList::validArgs.append("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "specify the source time"
    );

    argList::addOption
    (
        "sourceRegion",
        "word",
        "specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "specify the target region"
    );
//     argList::addBoolOption
//     (
//         "parallelSource",
//         "the source is decomposed"
//     );

    argList::addBoolOption
    (
        "extend",
        "Use the closest cell value if a given cell is outside of the mesh"
    );
    argList::addOption
    (
        "maxR",
        "scalar|'maximum radius'",
        "Cut off radius to map"
    );
    argList::addOption
    (
        "centre",
        "vector|'(0 0 0)'",
        "Location of center in the target mesh"
    );
    argList::addOption
    (
        "sourceCentre",
        "vector|'(0 0 0)'",
        "Location of center in the source mesh"
    );
    argList::addOption
    (
        "additionalFields",
        "wordList|'(rho U)'",
        "List of additional fields to map"
    );
    argList::addBoolOption
    (
        "uniform",
        "Copy uniform objects (not time)"
    );

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    //- Select time
    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);


    #include "createNamedMesh.H"

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
    word sourceRegion = fvMesh::defaultRegion;
    if (args.optionFound("sourceRegion"))
    {
        sourceRegion = args["sourceRegion"];
        Info<< "Source region: " << sourceRegion << endl;
    }

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = fvMesh::defaultRegion;
    if (args.optionFound("targetRegion"))
    {
        targetRegion = args["targetRegion"];
        Info<< "Target region: " << targetRegion << endl;
    }

//     const bool parallelSource = args.optionFound("parallelSource");

    scalar maxR(-1);
    if (args.optionFound("maxR"))
    {
        maxR = args.optionRead<scalar>("maxR");
    }
    else if (args.optionFound("extend"))
    {
        maxR = great;
    }


    wordList additionalFieldsNames;
    if (args.optionFound("additionalFields"))
    {
        additionalFieldsNames = args.optionRead<wordList>("additionalFields");
    }
    HashSet<word> additionalFields(additionalFieldsNames);

    bool copyUniform = args.optionFound("uniform");

    fvMesh& targetMesh = mesh;
    Time& runTimeTarget = runTime;

    const string caseDirOrig = getEnv("FOAM_CASE");
    const string caseNameOrig = getEnv("FOAM_CASENAME");
    setEnv("FOAM_CASE", rootDirSource/caseDirSource, true);
    setEnv("FOAM_CASENAME", caseDirSource, true);
    Time runTimeSource
    (
        Time::controlDictName,
        rootDirSource,
        caseDirSource
    );
    setEnv("FOAM_CASE", caseDirOrig, true);
    setEnv("FOAM_CASENAME", caseNameOrig, true);

//     if (parallelSource)
//     {
//         IOdictionary decompositionDict
//         (
//             IOobject
//             (
//                 "decomposeParDict",
//                 runTimeSource.system(),
//                 runTimeSource,
//                 IOobject::MUST_READ_IF_MODIFIED,
//                 IOobject::NO_WRITE
//             )
//         );
//         label nProcs = readInt(decompositionDict.lookup("numberOfSubdomains"));
//
//         Info<< "Create target mesh\n" << endl;
//         Info<< "Target mesh size: " << meshTarget.nCells() << endl;
//
//         for (int proci=0; proci<nProcs; proci++)
//         {
//             Info<< nl << "Source processor " << proci << endl;
//             fileName parCaseDirSource
//             (
//                 caseDirSource/fileName(word("processor") + name(proci))
//             );
//             Info<<rootDirSource<<endl;
//
//             Time runTimeSource
//             (
//                 Time::controlDictName,
//                 rootDirSource,
//                 caseDirSource/fileName(word("processor") + name(proci))
//             );
//             Pout<<"time"<<endl;
//
//             #include "setTimeIndex.H"
//
//             fvMesh meshSource
//             (
//                 IOobject
//                 (
//                     sourceRegion,
//                     runTimeSource.timeName(),
//                     runTimeSource
//                 )
//             );
//             Info<< "Source mesh size: " << meshSource.nCells() << tab
//                 << "Target mesh size: " << meshTarget.nCells() << nl << endl;
//
//             Vector<label> sourceGeoD(meshSource.geometricD());
//             Vector<label> sourceSolD(meshSource.solutionD());
//             Vector<label> sourceD(0, 0, 0);
//             forAll(sourceGeoD, i)
//             {
//                 if (sourceGeoD[i] == sourceSolD[i])
//                 {
//                     sourceD[i] = 1;
//                 }
//             }
//
//             vector rotDir(0.0, 0.0, 0.0);
//             forAll(sourceD, i)
//             {
//                 if (sourceD[i] != targetD[i])
//                 {
//                     rotDir[i] = 1.0;
//                 }
//             }
//             vector solDir(0.0, 0.0, 0.0);
//             forAll(sourceD, i)
//             {
//                 if (sourceD[i] == 1 && targetD[i] == 1)
//                 {
//                     solDir[i] = 1.0;
//                 }
//             }
//
//             vector centre(sum(cmptMultiply(meshSource.cellCentres(), rotDir)));
//             centre = cmptMultiply(centre, rotDir);
//
//             // Create map
//             labelList cellMap(meshTarget.nCells(), -1);
//             vectorField n(meshTarget.nCells(), Zero);
//             forAll(cellMap, celli)
//             {
//                 vector ptTarget =
//                     cmptMultiply(meshTarget.C()[celli], solDir)
//                   + cmptMultiply(meshTarget.C()[celli], rotDir);
//                 scalar r = mag(ptTarget);
//                 vector ptSource(r*solDir + centre);
//                 cellMap[celli] = meshSource.findCell(ptSource);
//                 n[celli] =
//                     ((ptSource - centre) & (solDir + rotDir))
//                    *(solDir + rotDir);
//             }
//
//
//             mapFields(meshSource, meshTarget, cellMap, n);
//         }
//     }
//     else
    {
        Time runTimeSource
        (
            Time::controlDictName,
            rootDirSource,
            caseDirSource
        );
        #include "setTimeIndex.H"

        Info<< "Create meshes\n" << endl;

        fvMesh sourceMesh
        (
            IOobject
            (
                sourceRegion,
                runTimeSource.timeName(),
                runTimeSource
            )
        );

        vector sourceCentre(Zero);
        if (args.optionFound("sourceCentre"))
        {
            sourceCentre = args.optionRead<vector>("sourceCentre");
        }
        else
        {
            sourceCentre = calculateCentre(sourceMesh);
        }

        vector targetCentre(sourceCentre);
        if (args.optionFound("centre"))
        {
            targetCentre = args.optionRead<vector>("centre");
        }

        Pair<vector> sourceAxis(calculateAxis(sourceMesh));
        Pair<vector> targetAxis(calculateAxis(targetMesh));
        vector rotationAxis = sourceAxis[1] - targetAxis[1];
        vector rAxis = sourceAxis[0];

        Info<< "Source mesh size: " << sourceMesh.nCells() << tab
            << "Target mesh size: " << targetMesh.nCells() << nl << endl;

        labelList cellMap;
        tensorField R;
        calcMapAndR
        (
            sourceMesh,
            targetMesh,
            maxR,
            targetCentre,
            sourceCentre,
            rotationAxis,
            rAxis,
            cellMap,
            R
        );

        // Map fields from the source mesh to the target mesh
        mapFields
        (
            sourceMesh,
            targetMesh,
            cellMap,
            R,
            additionalFieldsNames
        );

        if (copyUniform)
        {
            Info<<"here"<<endl;
            IOobjectList uniformObjects
            (
                sourceMesh,
                sourceMesh.time().timeName()/"uniform"
            );
            forAllConstIter
            (
                IOobjectList,
                uniformObjects,
                iter
            )
            {
                Info<<iter()->name()<<endl;
                Info<<iter()->objectPath()<<endl;
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
