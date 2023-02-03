/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
    Adds and empty patch to mesh and optionally fields

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "dynMeshTools.H"
#include "processorPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patch name");
    argList::validArgs.append("patch type");


    //- Add options
    timeSelector::addOptions(true, false);

    argList::addBoolOption
    (
        "noFields",
        "Do not add internal patch to fields"
    );
    argList::addBoolOption
    (
        "overwrite",
        "Write the mesh to constant"
    );
    argList::addBoolOption
    (
        "first",
        "add the new patch as the first patch"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Select time
    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    // Store face instance
    const word oldFacesInstance = mesh.facesInstance();

    // Read arguments
    const word patchName(args[1]);
    const word patchType(args[2]);

    // Do not write fields
    // Usefull if a refined mesh is needed before mesh manipulation
    bool noFields(args.optionFound("noFields"));
    bool overwrite(args.optionFound("overwrite"));

    if (!noFields)
    {
        meshTools::readAndStoreFields(mesh);
    }

    // Previous boundary patches
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find patch ID of specified patch
    label patchID = pbm.findPatchID(patchName);
    label startFace = mesh.nInternalFaces();

    if (patchID != -1)
    {
        WarningInFunction
            << patchName << " already exists!" << endl;
        return 0;
    }


    if (args.optionFound("first"))
    {
        patchID = 0;
    }
    else if (Pstream::parRun())
    {
        forAll(pbm, patchi)
        {
            startFace = pbm[patchi].start();
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                patchID = patchi;
                break;
            }
        }
    }
    else
    {
        patchID = pbm.size();
    }

    autoPtr<polyPatch> newPatch
    (
        polyPatch::New
        (
            patchType,
            patchName,
            0,
            startFace,
            patchID,
            pbm
        )
    );
    mesh.addPatch
    (
        patchID,
        newPatch(),
        dictionary(),
        patchType,
        true
    );
    Info<< "Added " << patchName << " at index " << patchID << nl << endl;


    // Write mesh and cell levels
    if (overwrite)
    {
        mesh.setInstance(oldFacesInstance);
    }
    else
    {
        runTime++;
        mesh.setInstance(runTime.timeName());
    }

    Info<<"Writing mesh to " << mesh.facesInstance() << nl << endl;
    mesh.write();

    Info<< "\nEnd\n" << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
