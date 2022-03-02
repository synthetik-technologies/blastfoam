/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Calculate distance from 0 levelSet using searchableSurfaces, and set the
    associated volumeFractionField. Optional refinement of interface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "backupSearchableSurface.H"
#include "systemDict.H"
#include "calcAngleFraction.H"
#include "dynMeshTools.H"

#include "fvMeshHexRefiner.H"
#include "topoSetList.H"
#include "levelSetModel.H"
#include "IOobjectList.H"
#include "fvc.H"

using namespace Foam;

void setPhase
(
    autoPtr<fvMeshRefiner>& refiner,
    fvMesh& mesh,
    volScalarField& alpha,
    levelSetModel& LSModel,
    const dictionary& dict
)
{
    const scalar angleFraction = calcAngleFraction(mesh);

    //- List of sources (and backups if present)
    // Stored to reduce the number of reads
    PtrList<backupSearchableSurface> regions
    (
        dict.lookup("initialRegions"),
        backupSearchableSurface::iNew(mesh)
    );

    // Collection of searchable surfaces
    labelList levels(regions.size());
    UPtrList<searchableSurface> surfaces(regions.size());
    forAll(regions, regionI)
    {
        levels[regionI] =
            regions[regionI].dict().lookupOrDefault("level", 0);
        surfaces.set
        (
            regionI,
            &const_cast<searchableSurface&>(regions[regionI].surface())
        );
    }

    label maxLevel = max(levels);

    // Maximum number of iterations
    label iter = 0;
    label maxIter =
        max
        (
            1,
            max
            (
                2*maxLevel,
                gMax(refiner->cellLevel())*2
            )
        );

    // Flag for final iteration
    bool end = false;

    // Flag to initiate end
    bool prepareToStop = (maxIter == 1);

    volScalarField error
    (
        IOobject
        (
            "error",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        -1.0
    );

    while(!end)
    {
        if (maxIter <= iter)
        {
            prepareToStop = true;
            WarningInFunction << "Reach maximum number of iteration" << endl;
        }

        error == -1.0;

        // Check if this is the final iteration so correct cell shapes are used
        if (prepareToStop)
        {
            end = true;
        }

        // Read default values and set fields
        alpha == dict.lookupOrDefault("defaultValue", 0.0);

        // List of saved cells (per cell set)
        labelListList savedCells(regions.size());
        labelListList savedPoints(regions.size());

        Info<< "Setting field region values" << endl;
        forAll(regions, regionI)
        {
            const dictionary& regionDict = regions[regionI].dict();
            regions[regionI].allowBackup(!end);
            savedCells[regionI] = regions[regionI].selectPoints
            (
                mesh.cellCentres()
            );
            const labelList& selectedCells = savedCells[regionI];

            if (!returnReduce(selectedCells.size(), sumOp<label>()))
            {
                WarningInFunction
                    << "No cells were selected for using " << nl
                    << regionDict
                    << "To expand searchable region add backup " << nl
                    << "Region or expand backup region." << endl;
            }

            bool set = selectedCells.size();
            reduce(set, orOp<bool>());

            if (set)
            {
                Info<< "    Selected "
                    << returnReduce(selectedCells.size(), sumOp<label>())
                    << " cells" << endl;

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
                Info<< nl << endl;

                bool add = !regionDict.lookupOrDefault("remove", false);
                if (add)
                {
                    UIndirectList<scalar>(alpha, selectedCells) = 1.0;
                }
                else
                {
                    UIndirectList<scalar>(alpha, selectedCells) = 0.0;
                }
            }
        }

        // Update boundary conditions of fields used for refinement
        alpha.correctBoundaryConditions();

        // Update error and mesh if not the final iteration
        if (refiner.valid())
        {
            LSModel.levelSet() = LSModel.calcLevelSet(alpha, surfaces);
            alpha = LSModel.alpha();

            forAll(error, celli)
            {
                if (alpha[celli] < 0.99 && alpha[celli] > 0.01)
                {
                    error[celli] = 1.0;
                }
            }

            labelList maxCellLevel(mesh.nCells(), -1);
            forAll(regions, regionI)
            {
                // Set specified cells to be refined
                labelList refineCells;
                labelList refinePoints;
                if (dict.lookupOrDefault("refineInternal", false))
                {
                    refineCells = savedCells[regionI];
                }
                else
                {
                    refineCells =
                        topoSetList::extractInterfaceCells
                        (
                            mesh,
                            savedCells[regionI]
                        );
                }

                if (dict.lookupOrDefault("refinePoints", false))
                {
                    labelList selectedPoints
                    (
                        regions[regionI].selectPoints(mesh.points())
                    );
                    refinePoints = topoSetList::extractSelectedPoints
                    (
                        mesh,
                        dict,
                        selectedPoints,
                        false
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
                const labelListList& pointCells = mesh.pointCells();
                if (overwriteLevel)
                {
                    forAll(refineCells, celli)
                    {
                        maxCellLevel[refineCells[celli]] = levels[regionI];
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
            // Update mesh (return if mesh changes)
            if (!end)
            {
                prepareToStop = !refiner->refine(error, maxCellLevel);
            }
        }
        iter++;
    }

    LSModel.levelSet() = LSModel.calcLevelSet(alpha, surfaces);
    alpha = LSModel.alpha();
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

    // Store original mesh instance
    const fileName oldFacesInstance = mesh.facesInstance();

    const dictionary levelSetProperties
    (
        systemDict("levelSetProperties", args, mesh)
    );

    // Update All fields in the current time folder
    // Resizes to make sure all fields are consistent with the mesh
    bool updateAll(args.optionFound("updateAll"));
    bool noRefine(args.optionFound("noRefine"));
    bool overwrite(args.optionFound("overwrite"));
    bool noHistory(args.optionFound("noHistory"));

    //- Is the mesh balanced
    autoPtr<fvMeshRefiner> refiner;
    if (!noRefine)
    {
        dictionary refineDict
        (
            levelSetProperties.optionalSubDict
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

    // Collect the phases to set
    wordList phases(levelSetProperties.lookup("phases"));

    PtrList<volScalarField> alphas(phases.size());
    PtrList<levelSetModel> LSModels(phases.size());
    label phaseI = 0;
    forAll(phases, phasei)
    {
        IOobject alphaIO
        (
            IOobject::groupName("alpha", phases[phasei]),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if (alphaIO.typeHeaderOk<volScalarField>(true))
        {
            alphas.set
            (
                phaseI,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("alpha", phases[phasei]),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
            LSModels.set
            (
                phaseI,
                new levelSetModel
                (
                    alphas[phasei],
                    levelSetProperties.subDict(phases[phasei]),
                    false
                )
            );
            phases[phaseI] = phases[phasei];
            phaseI++;
        }
    }
    alphas.resize(phaseI);
    LSModels.resize(phaseI);
    phases.resize(phaseI);

    // Read in all fields to allow resizing
    if (updateAll)
    {
        meshTools::readAndStoreFields(mesh);
    }

    //- List of sources (and backups if present)
    // Stored to reduce the number of reads
    forAll(phases, phasei)
    {
        Info<<"Creating levelSet function for " << phases[phasei] << endl;
        setPhase
        (
            refiner,
            mesh,
            alphas[phasei],
            LSModels[phasei],
            levelSetProperties.subDict(phases[phasei])
        );
        Info<< endl;
    }

    if (overwrite)
    {
        mesh.setInstance(oldFacesInstance);
    }

    bool writeMesh = false;
    if (refiner.valid())
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

    if (noHistory)
    {
        refiner.clear();
    }

    // Write all fields
    runTime.write();
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
