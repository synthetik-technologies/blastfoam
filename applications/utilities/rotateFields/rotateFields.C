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
#include "labelVector.H"
#include "wedgeFvPatch.H"
#include "IOobjectList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void mapVolFields
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const labelList& cellMap,
    const IOobjectList& objects,
    const tensorField& R
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    IOobjectList fields = objects.lookupClass(fieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            Info<< "    interpolating " << fieldIter()->name() << endl;

            // Read field fieldSource
            fieldType fieldSource(*fieldIter(), meshSource);

            // Read fieldTarget
            fieldType fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

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

template<>
void mapVolFields<scalar>
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const labelList& cellMap,
    const IOobjectList& objects,
    const tensorField&
)
{
    typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
    IOobjectList fields = objects.lookupClass(fieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            Info<< "    interpolating " << fieldIter()->name() << endl;

            // Read field fieldSource
            fieldType fieldSource(*fieldIter(), meshSource);

            // Read fieldTarget
            fieldType fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

            forAll(cellMap, celli)
            {
                label cellj = cellMap[celli];
                if (cellj != -1)
                {
                    fieldTarget[celli] = fieldSource[cellj];
                }
            }
            fieldTarget.write();
        }
    }
}

void mapFields
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const labelList& cellMap,
    const HashSet<word>& fields,
    const tensorField& R
)
{
    IOobjectList objects(meshSource, meshSource.time().timeName());
    if (fields.toc().size() > 0)
    {
        forAllIter(IOobjectList, objects, objectIter)
        {
            if (!fields.found(objectIter()->name()))
            {
                objects.erase(objectIter);
            }
        }
    }


    mapVolFields<scalar>(meshSource, meshTarget, cellMap, objects, tensorField());
    mapVolFields<vector>(meshSource, meshTarget, cellMap, objects, R);
//     mapVolFields<sphericalTensor>(meshSource, meshTarget, cellMap, objects, R);
//     mapVolFields<symmTensor>(meshSource, meshTarget, cellMap, objects, R);
    mapVolFields<tensor>(meshSource, meshTarget, cellMap, objects, R);
}


void calcMapAndR
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const scalar& maxR,
    const vector& centreTarget,
    labelList& cellMap,
    tensorField& R
)
{
    Vector<label> targetGeoD(meshTarget.geometricD());
    Vector<label> targetSolD(meshTarget.solutionD());
    Vector<label> targetD(0, 0, 0);
    forAll(targetGeoD, i)
    {
        if (targetGeoD[i] == targetSolD[i])
        {
            targetD[i] = 1;
        }
    }

    Vector<label> sourceGeoD(meshSource.geometricD());
    Vector<label> sourceSolD(meshSource.solutionD());
    Vector<label> sourceD(0, 0, 0);
    forAll(sourceGeoD, i)
    {
        if (sourceGeoD[i] == sourceSolD[i])
        {
            sourceD[i] = 1;
        }
    }

    vector rotDir(0.0, 0.0, 0.0);
    forAll(sourceD, i)
    {
        if (sourceD[i] != targetD[i])
        {
            rotDir[i] = 1.0;
        }
    }
    vector solDir(0.0, 0.0, 0.0);
    forAll(sourceD, i)
    {
        if (sourceD[i] == 1 && targetD[i] == 1)
        {
            solDir[i] = 1.0;
        }
    }

    vector centre
    (
        (
            sum
            (
                cmptMultiply
                (
                    meshSource.C()()*meshSource.V(),
                    rotDir
                )
            )/sum(meshSource.V())
        ).value()
    );

    // Create map
    cellMap = labelList(meshTarget.nCells(), -1);
    R = tensorField(meshTarget.nCells(), tensor::I);
    forAll(cellMap, celli)
    {
        // Get the position on the target mesh
        vector ptTarget =
            cmptMultiply(meshTarget.cellCentres()[celli], (solDir + rotDir));

        // Offset from the center of the source mesh
        vector nTarget = ptTarget - centre - centreTarget;

        // Radius
        scalar r = mag(nTarget);

        // Remove points outside of the maximum radius
        if (r > maxR)
        {
            r = great;
        }

        // Offset from the source mesh center (only solved directions)
        vector nSource(r*solDir);

        // Actual point on the source mesh
        vector ptSource = nSource + centre;

        // Map from the source mesh to the target mesh
        cellMap[celli] = meshSource.findCell(ptSource);

        // Extend radius is the target point is outside of the source mesh
        if (maxR < 0 && cellMap[celli] < 0)
        {
            cellMap[celli] = meshSource.findNearestCell(ptSource);
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
        "subtract",
        "subtract mapped source from target"
    );
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
        "fields",
        "wordList|'(rho U)'",
        "List of fields to map"
    );

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
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

    scalar maxR(great);
    if (args.optionFound("maxR"))
    {
        maxR = args.optionRead<scalar>("maxR");
    }
    else if (args.optionFound("extend"))
    {
        maxR = -great;
    }

    vector centreTarget(Zero);
    if (args.optionFound("centre"))
    {
        centreTarget = args.optionRead<vector>("centre");
    }

    wordList fieldNames;
    if (args.optionFound("fields"))
    {
        fieldNames = args.optionRead<wordList>("fields");
    }
    HashSet<word> fields(fieldNames);

    if (args.optionFound("parallel"))
    {
        caseDirTarget =
            caseDirTarget/fileName(word("processor")+ name(Pstream::myProcNo()));
    }

    fvMesh& meshTarget = mesh;
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

        fvMesh meshSource
        (
            IOobject
            (
                sourceRegion,
                runTimeSource.timeName(),
                runTimeSource
            )
        );

        Info<< "Source mesh size: " << meshSource.nCells() << tab
            << "Target mesh size: " << meshTarget.nCells() << nl << endl;

        labelList cellMap;
        tensorField R;
        calcMapAndR(meshSource, meshTarget, maxR, centreTarget, cellMap, R);

        // Map fields from the source mesh to the target mesh
        mapFields(meshSource, meshTarget, cellMap, fieldNames, R);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
