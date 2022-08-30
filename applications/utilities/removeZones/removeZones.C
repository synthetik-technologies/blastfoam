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
    Simple application to remove zones. Useful for removing the zones created
    by gmshToFoam

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "cellZone.H"
#include "faceZone.H"
#include "pointZone.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //- Add options
    timeSelector::addOptions(true, false);

    argList::addBoolOption
    (
        "clear",
        "Clear all cell zones"
    );
    argList::addBoolOption
    (
        "allCells",
        "Clear cellZones"
    );
    argList::addBoolOption
    (
        "allFaces",
        "Clear faceZones"
    );
    argList::addOption
    (
        "allPoints",
        "Clear pointZones"
    );
    argList::addOption
    (
        "cells",
        "List of cellZones to remove"
    );
    argList::addOption
    (
        "faces",
        "List of faceZones to remove"
    );
    argList::addOption
    (
        "points",
        "List of pointZones to remove"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Select time
    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const fileName meshInstance(mesh.facesInstance());

    //- Full clear of selected
    bool clear(args.optionFound("clear"));
    bool allCells(args.optionFound("allCells") || clear);
    bool allFaces(args.optionFound("allFaces") || clear);
    bool allPoints(args.optionFound("allPoints") || clear);

    List<cellZone*> newCellZones(mesh.cellZones().size());
    List<faceZone*> newFaceZones(mesh.faceZones().size());
    List<pointZone*> newPointZones(mesh.pointZones().size());

    label czi = 0;
    label fzi = 0;
    label pzi = 0;

    if (!allCells)
    {
        wordHashSet cellZonesToRemove;
        args.optionReadIfPresent("cells", cellZonesToRemove);

        const meshCellZones& cellZones = mesh.cellZones();
        forAll(cellZones, zonei)
        {
            const word& zoneName = cellZones[zonei].name();
            if (!cellZonesToRemove.found(zoneName))
            {
                newCellZones[czi] =
                    new cellZone
                    (
                        zoneName,
                        cellZones[zonei],
                        czi,
                        mesh.cellZones()
                    );
                czi++;
            }
        }
    }

    if (!allFaces)
    {
        wordHashSet faceZonesToRemove;
        args.optionReadIfPresent("faces", faceZonesToRemove);

        const meshFaceZones& faceZones = mesh.faceZones();
        forAll(faceZones, zonei)
        {
            const word& zoneName = faceZones[zonei].name();
            if (!faceZonesToRemove.found(zoneName))
            {
                newFaceZones[fzi] =
                    new faceZone
                    (
                        zoneName,
                        faceZones[zonei],
                        faceZones[zonei].flipMap(),
                        fzi,
                        mesh.faceZones()
                    );
                fzi++;
            }
        }
    }

    if (!allPoints)
    {
        wordHashSet pointZonesToRemove;
        args.optionReadIfPresent("points", pointZonesToRemove);

        const meshPointZones& pointZones = mesh.pointZones();
        forAll(pointZones, zonei)
        {
            const word& zoneName = pointZones[zonei].name();
            if (!pointZonesToRemove.found(zoneName))
            {
                newPointZones[pzi] =
                    new pointZone
                    (
                        zoneName,
                        pointZones[zonei],
                        pzi,
                        mesh.pointZones()
                    );
                pzi++;
            }
        }
    }

    newCellZones.setSize(czi);
    newFaceZones.setSize(fzi);
    newPointZones.setSize(pzi);

    mesh.pointZones().clear();
    mesh.faceZones().clear();
    mesh.cellZones().clear();

    mesh.addZones(newPointZones, newFaceZones, newCellZones);

    Info<< "writing mesh to " << meshInstance << endl;
    mesh.setInstance(meshInstance);
    mesh.write();

    Info<< "\nEnd\n" << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
