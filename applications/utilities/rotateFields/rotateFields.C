/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | SYnthetik Applied Technologies
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
    Rotate fields from 1-D to 2-D or 2-D to 3-D. Only for
    axisymmetric cases. Currently mapping is not conservative so only the
    nearest cell value (based on cell centres) will be used. Optionally,
    refinement can be used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "labelVector.H"
#include "wedgeFvPatch.H"
#include "IOobjectList.H"
#include "HashSet.H"
#include "UautoPtr.H"
#include "genericFvPatchField.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

#include "fvMeshRefiner.H"
#include "errorEstimator.H"

#include "mappingFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapFields
(
    const PtrList<fvMesh>& sourceMeshes,
    const fvMesh& targetMesh,
    const List<cellInfoList>& cellMap,
    const List<cellInfoList>& extendedCellMap,
    const tensorField& R,
    const HashSet<word>& additionalFields,
    const bool store
)
{
    Info<< "Mapping fields" << endl;
    IOobjectList objects(sourceMeshes[0], sourceMeshes[0].time().timeName());

    mapVolFields<scalar>
    (
        sourceMeshes,
        targetMesh,
        cellMap,
        extendedCellMap,
        objects,
        R,
        additionalFields,
        store
    );
    mapVolFields<vector>
    (
        sourceMeshes,
        targetMesh,
        cellMap,
        extendedCellMap,
        objects,
        R,
        additionalFields,
        store
    );
    mapVolFields<sphericalTensor>
    (
        sourceMeshes,
        targetMesh,
        cellMap,
        extendedCellMap,
        objects,
        R,
        additionalFields,
        store
    );
    mapVolFields<symmTensor>
    (
        sourceMeshes,
        targetMesh,
        cellMap,
        extendedCellMap,
        objects,
        R,
        additionalFields,
        store
    );
    mapVolFields<tensor>
    (
        sourceMeshes,
        targetMesh,
        cellMap,
        extendedCellMap,
        objects,
        R,
        additionalFields,
        store
    );
    Info<< endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Add options
    addOptions();

    #include "setRootCase.H"
    parRun = Pstream::parRun();

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();

    const bool nearest = args.optionFound("nearest");
    if (!isDir(casePath))
    {
        FatalErrorInFunction
            << casePath << " is not a valid directory" << endl
            << abort(FatalError);
    }

    Info<< "Source: " << casePath << " " << caseDirSource << endl;
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

    const bool parallelSource = args.optionFound("parallelSource");

    scalar maxR(-1);
    if (args.optionFound("maxR"))
    {
        maxR = args.optionRead<scalar>("maxR");
        Info<< "Maximum distance from target centre is " << maxR << endl;
    }
    else if (args.optionFound("extend"))
    {
        maxR = great;
        Info<< "Extending mapping to the edge of the domain" << endl;
    }

    wordList additionalFieldNames;
    if (args.optionFound("additionalFields"))
    {
        additionalFieldNames =
            args.optionRead<wordList>("additionalFields");
    }
    HashSet<word> additionalFields(additionalFieldNames);

    bool copyUniform = args.optionFound("uniform");


    const string caseDirOrig = getEnv("FOAM_CASE");
    const string caseNameOrig = getEnv("FOAM_CASENAME");

    Time targetRunTime(Foam::Time::controlDictName, args);
    fvMesh targetMesh
    (
        IOobject
        (
            targetRegion,
            targetRunTime.timeName(),
            targetRunTime,
            IOobject::MUST_READ
        )
    );
    Info<<"Created target mesh"<<endl;

    const polyMesh::cellDecomposition decompMode =
        args.optionFound("tets") ? polyMesh::CELL_TETS : polyMesh::FACE_DIAG_TRIS;

    PtrList<Time> sourceRunTimes;
    PtrList<fvMesh> sourceMeshes;
    label nSourceCells = 0;

    // Create source case argList
    argList sourceArgs(args);
    const_cast<ParRunControl&>(sourceArgs.parRunControl()) = ParRunControl();

    if (parallelSource)
    {
        label nProcs = fileHandler().nProcs(rootDirSource/caseDirSource);
        reduce(nProcs, maxOp<label>());
        if (nProcs < 1)
        {
            FatalErrorInFunction
                << "Trying to map from a parallel case, but no processor" << nl
                << "directories were found. remove the \"parallelSource\"" << nl
                << "for serial cases." << endl
                << abort(FatalError);
        }

        setParRun(false);
        sourceRunTimes.setSize(nProcs);
        sourceMeshes.setSize(nProcs);

        setEnv("FOAM_CASE", rootDirSource/caseDirSource, true);
        setEnv("FOAM_CASENAME", caseDirSource, true);
        for (int proci=0; proci < nProcs; proci++)
        {
            sourceRunTimes.set
            (
                proci,
                new Time
                (
                    rootDirSource,
                    caseDirSource/fileName(word("processor") + name(proci))
                )
            );
            Time& runTimeSource = sourceRunTimes[proci];
            const_cast<dictionary&>(runTimeSource.controlDict()) =
                targetRunTime.controlDict();
            #include "setTimeIndex.H"

            sourceMeshes.set
            (
                proci,
                new fvMesh
                (
                    IOobject
                    (
                        sourceRegion,
                        runTimeSource.timeName(),
                        runTimeSource
                    )
                )
            );
            nSourceCells += sourceMeshes[proci].nCells();
        }
        resetParRun();

        Info<< nl << "Read " << nProcs << " source processor meshes" << endl;
    }
    else
    {
        sourceRunTimes.setSize(1);
        sourceMeshes.setSize(1);

        sourceRunTimes.set
        (
            0,
            new Time
            (
                Time::controlDictName,
                rootDirSource,
                caseDirSource
            )
        );
        Time& runTimeSource = sourceRunTimes[0];
        #include "setTimeIndex.H"

        sourceMeshes.set
        (
            0,
            new fvMesh
            (
                IOobject
                (
                    sourceRegion,
                    runTimeSource.timeName(),
                    runTimeSource
                )
            )
        );
        Info<< "Created source mesh\n" << endl;
        nSourceCells += sourceMeshes[0].nCells();
    }
    setEnv("FOAM_CASE", caseDirOrig, true);
    setEnv("FOAM_CASENAME", caseNameOrig, true);

    Info<< "\nSource time: " << sourceRunTimes[0].value()
        << "\nTarget time: " << targetRunTime.value()
        << nl << endl;

    vector sourceSumCV = Zero;
    scalar sourceSumV = 0.0;
    PtrList<indexedOctree<treeDataCell>> icos(sourceMeshes.size());
    setParRun(false);

    forAll(sourceMeshes, proci)
    {
        const fvMesh& sourceMesh = sourceMeshes[proci];
        sourceSumCV += sum(sourceMesh.C()*sourceMesh.V()).value();
        sourceSumV += sum(sourceMesh.V()).value();

        treeBoundBox meshBb(sourceMesh.bounds());

        // Calculate typical cell related size to shift bb by.
        scalar typDim = meshBb.avgDim()/(2.0*Foam::cbrt(scalar(sourceMesh.nCells())));

        treeBoundBox shiftedBb
        (
            meshBb.min(),
            meshBb.max() + vector(typDim, typDim, typDim)
        );

        icos.set
        (
            proci,
            new indexedOctree<treeDataCell>
            (
                treeDataCell(true, sourceMesh, decompMode),
                shiftedBb,
                10,         // maxLevel
                100,        // leafsize
                10.0        // duplicity
            )
        );
        sourceMesh.tetBasePtIs();
    }
    resetParRun();
    Info<< "created source meshes" << nl << endl;

    Pair<vector> sourceAxis(calculateAxis(sourceMeshes[0]));
    Pair<vector> targetAxis(calculateAxis(targetMesh));
    vector rotationAxis = sourceAxis[1] - targetAxis[1];
    vector rAxis = sourceAxis[0];

    vector sourceCentre = cmptMultiply(sourceSumCV, sourceAxis[1])/sourceSumV;
    vector targetCentre(sourceCentre);
    if (args.optionFound("centre"))
    {
        targetCentre = args.optionRead<vector>("centre");
    }

    Info<< "Source centre: " << sourceCentre << nl
        << "Target centre: " << targetCentre << endl;

    Info<< "Source mesh size: " << nSourceCells << endl;

    if (!args.optionFound("refine"))
    {
        Info<< "Target mesh size: " << targetMesh.nCells() << nl << endl;
        List<cellInfoList> cellMap(targetMesh.nCells());
        List<cellInfoList> extendedCellMap(targetMesh.nCells());
        tensorField R(targetMesh.nCells(), tensor::I);

        calcMapAndR
        (
            icos,
            sourceMeshes,
            targetMesh,
            maxR,
            sourceCentre,
            targetCentre,
            rotationAxis,
            rAxis,
            nearest,
            cellMap,
            extendedCellMap,
            R
        );

        // Map fields from the source mesh to the target mesh
        mapFields
        (
            sourceMeshes,
            targetMesh,
            cellMap,
            extendedCellMap,
            R,
            additionalFieldNames
        );
    }
    else
    {
        refine
        (
            icos,
            sourceMeshes,
            targetMesh,
            maxR,
            sourceCentre,
            targetCentre,
            rotationAxis,
            rAxis,
            nearest,
            additionalFieldNames
        );
        targetRunTime.writeNow();
    }

    if (copyUniform)
    {
        fileName local = "uniform";
        fileName path = targetMesh.time().timePath();

        IOobjectList uniformObjects
        (
            sourceMeshes[0],
            sourceRunTimes[0].timeName()/local
        );
        forAllConstIter
        (
            IOobjectList,
            uniformObjects,
            iter
        )
        {
            fileName name = iter()->name();
            if (name != "time")
            {
                fileName srcPath = iter()->objectPath();
                cp
                (
                    iter()->objectPath(),
                    path/local/name
                );
            }
        }
    }

    // Write points0 field to time directory
    if (args.optionFound("points0"))
    {
        Info<< "Writing points0" << endl;
        pointIOField
        (
            IOobject
            (
                "points0",
                targetMesh.facesInstance(),
                polyMesh::meshSubDir,
                targetMesh
            ),
            targetMesh.points()
        ).write();
    }

    Info<< nl << "Finished" << endl
        << "ExecutionTime = " << targetRunTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << targetRunTime.elapsedClockTime() << " s" << nl
        << endl;

    return 0;
}


// ************************************************************************* //
