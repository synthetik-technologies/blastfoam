/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
05-08-2020 Jeff Heylmun     | Made blastFSIFoam solver
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

Application
    blastFSIFoam

Description

    EXPERIMENTAL

    Solver for solving couple fluid structure interaction using the
    standard blastFoam classes and the solidDisplacementFoam solver.
    The fluid phase uses moving mesh, but the solid phase does not. For
    this reason modified mapped boundaries use an offset of the solid
    displacement to correctly map the boundary information between regions.

    This solver is currently under developments and is not stable.
    Adaptive refinement does work, however it is not robust. Dynamic
    load balancing should not be used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "timeIntegrator.H"
#include "phaseCompressibleSystem.H"
#include "compressibleCourantNo.H"
#include "fluidThermoModel.H"
#include "solidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "blastCompressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "fvOptions.H"

#include "mappedPatchSelector.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "setBoundaryDisplacementFields.H"

    while (runTime.run())
    {
        #include "refineMeshes.H"

        #include "readTimeControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "updateMeshes.H"

        #include "clearPatches.H"

        // Solve
        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "solveFluid.H"
        }

        forAll(solidRegions, i)
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"
            #include "readSolidTimeControls.H"

            #include "solveSolid.H"
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
