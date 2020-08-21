/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
20-06-2020 Jeff Heylmun     | Added use of blastFoam thermodynamics
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
    blastMultiRegionFoam

Description
    Solver for transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, buoyancy effects, turbulence,
    and radiation modeling. Riemann fluxes are used to transport the fluid
    phase.

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
#include "mappedWallFvPatch.H"
#include "mappedPatchBase.H"


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

    while (runTime.run())
    {
        #include "readTimeControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //- Update internal error fields of all regions
        forAll(fluidRegions, i)
        {
            fluidRegions[i].updateError();
        }
        forAll(solidRegions, i)
        {
            solidRegions[i].updateError();
        }

        //- Update error field boundaries
        //  Done after all region errors have been updated to make sure
        //  all fields are up to date
        forAll(fluidRegions, i)
        {
            fluidRegions[i].updateErrorBoundaries();
        }
        forAll(solidRegions, i)
        {
            solidRegions[i].updateErrorBoundaries();
        }

        //- Update Meshes and check if balancing has occurred
        boolList balanced(fluidRegions.size() + solidRegions.size(), false);
        bool needRemap = false;
        label j = 0;
        forAll(fluidRegions, i)
        {
            fluidRegions[i].update();

            // Already cleared mapped patches
            balanced[j++] = !fluidRegions[i].balanced();
            needRemap = needRemap || fluidRegions[i].balanced();
        }
        forAll(solidRegions, i)
        {
            solidRegions[i].update();

            // Already cleared mapped patches
            balanced[j++] = !solidRegions[i].balanced();
            needRemap = needRemap || solidRegions[i].balanced();
        }

        // Clear mapped boundaries if one region has been balanced
        // Balanced meshes have already had their maps cleared
        if (needRemap)
        {
            j = 0;
            forAll(fluidRegions, i)
            {
                if (!balanced[j++])
                {
                    forAll(fluidRegions[i].boundaryMesh(), patchi)
                    {
                        if
                        (
                            isA<mappedWallFvPatch>
                            (
                                fluidRegions[i].boundary()[patchi]
                            )
                        )
                        {
                            polyBoundaryMesh& pbMesh =
                                const_cast<polyBoundaryMesh&>
                                (
                                    fluidRegions[i].boundaryMesh()
                                );
                            refCast<mappedPatchBase>(pbMesh[patchi]).clearOut();
                        }
                    }
                }
            }
            forAll(solidRegions, i)
            {
                if (!balanced[j++])
                {
                    forAll(solidRegions[i].boundaryMesh(), patchi)
                    {
                        if
                        (
                            isA<mappedWallFvPatch>
                            (
                                solidRegions[i].boundary()[patchi]
                            )
                        )
                        {
                            polyBoundaryMesh& pbMesh =
                                const_cast<polyBoundaryMesh&>
                                (
                                    solidRegions[i].boundaryMesh()
                                );
                            refCast<mappedPatchBase>(pbMesh[patchi]).clearOut();
                        }
                    }
                }
            }
        }

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
