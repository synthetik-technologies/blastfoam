/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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
    Converts case files from OpenFOAM-9 format to OpenFOAM-12 format

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeIfNew(IOdictionary& dict)
{
    if (!dict.headerOk())
    {
        dict.regIOobject::write();
    }
    else
    {
        WarningInFunction
            << dict.objectPath(false)
            << " already exists" << endl;
    }
}

int main(int argc, char *argv[])
{
    argList::addOption("solver", "Name of solver");
    argList::addBoolOption("multiphase", "blastEulerFoam case");

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"

    const word solver(args.optionLookupOrDefault<word>("solver", "blastFoam"));
    bool isBlast = solver == "blastFoam";
    bool isBlastEuler = solver == "blastFoamEuler" || args.optionFound("multiphase");

    typeIOobject<IOdictionary> phasePropertiesIO
    (
        "phaseProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    );
    if (phasePropertiesIO.headerOk())
    {
        IOdictionary phaseProperties(phasePropertiesIO);
        IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                runTime
            ),
            phaseProperties
        );
        if (isBlastEuler)
        {
            const wordList phases(phaseProperties.lookup("phases"));
            forAll(phases, i)
            {
                IOdictionary phasePhysicalProperties
                (
                    IOobject
                    (
                        IOobject::groupName("physicalProperties", phases[i]),
                        runTime.constant(),
                        runTime
                    ),
                    phaseProperties.subDict(phases[i])
                );
                writeIfNew(phasePhysicalProperties);

                typeIOobject<IOdictionary> turbulencePropertiesIO
                (
                    IOobject
                    (
                        IOobject::groupName("turbulenceProperties", phases[i]),
                        runTime.constant(),
                        runTime,
                        IOobject::MUST_READ
                    )
                );
                if (turbulencePropertiesIO.headerOk())
                {
                    IOdictionary turbulenceProperties
                    (
                        turbulencePropertiesIO
                    );
                    turbulenceProperties.rename
                    (
                        IOobject::groupName("momentumTransport", phases[i])
                    );
                    writeIfNew(turbulenceProperties);
                }
            }
        }
        else
        {
            writeIfNew(physicalProperties);

            typeIOobject<IOdictionary> turbulencePropertiesIO
            (
                IOobject
                (
                    "turbulenceProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ
                )
            );
            if (turbulencePropertiesIO.headerOk())
            {
                IOdictionary turbulenceProperties
                (
                    turbulencePropertiesIO
                );
                turbulenceProperties.rename("momentumTransport");
                writeIfNew(turbulenceProperties);
            }
        }

    }

    typeIOobject<IOdictionary> dynamicMeshDictIO
    (
        "dynamicMeshDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    );
    if (dynamicMeshDictIO.headerOk())
    {
        IOdictionary dynamicMeshDict0(dynamicMeshDictIO);

        if (dynamicMeshDict0.found("dynamicFvMesh"))
        {
            dynamicMeshDict0.rename(dynamicMeshDictIO.name() + "_old");
            writeIfNew(dynamicMeshDict0);

            dynamicMeshDictIO.readOpt() = IOobject::NO_READ;
            IOdictionary dynamicMeshDict(dynamicMeshDictIO);

            const word type(dynamicMeshDict0.lookup("dynamicFvMesh"));
            dynamicMeshDict0.remove("dynamicFvMesh");

            bool balance = false;
            bool refine = false;
            bool moving = false;
            word refiner = "hexRefiner";
            if (type == "adaptiveFvMesh" || type == "movingAdaptiveFvMesh")
            {
                refine = true;
                dynamicMeshDict0.readIfPresent("refiner", refiner);

                balance = true;
                dynamicMeshDict0.readIfPresent("balance", balance);
            }
            if (type == "dynamicMotionSolverFvMesh" || type == "movingAdaptiveFvMesh")
            {
                moving = true;
            }

            if (moving)
            {
                dictionary moverDict(dynamicMeshDict0);
                moverDict.set("type", "motionSolver");
                moverDict.set
                (
                    "libs",
                    fileNameList
                    (
                        {
                            "libfvMeshMoversMotionSolver.so",
                            "libfvMotionSolvers.so"
                        }
                    )
                );
                dynamicMeshDict.set("mover", moverDict);
            }

            if (refine)
            {
                dictionary topoChangerDict(dynamicMeshDict0);
                topoChangerDict.set("type", refiner);
                topoChangerDict.set
                (
                    "libs",
                    fileNameList({"libblastFvMeshTopoChangers.so"})
                );
                dynamicMeshDict.set("topoChanger", topoChangerDict);
            }
            if (balance)
            {
                dictionary distributorDict(dynamicMeshDict0);
                distributorDict.set("type", "redistributor");
                distributorDict.set
                (
                    "libs",
                    fileNameList({"libblastFvMeshDistributors.so"})
                );
                dynamicMeshDict.set("distributor", distributorDict);
            }
            dynamicMeshDict.regIOobject::write();
        }
    }

}


// ************************************************************************* //
