/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2025
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

Application
    blastEulerFoam

Description
    Multiphase compressible solver that uses Riemann solver to construct
    hyperbolic fluxes. Phases have unique velocities, internal energies
    and pressures.


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "phaseSystem.H"
#include "wedgeFvPatch.H"
#include "fvTimeIntegrator.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    scalar CoNum = 0.0;
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    maxCo = min(maxCo, integrator.maxCo());
    scalar mDotCoNum = 0.0;
    #include "EigenCourantNos.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        integrator.preUpdateMesh();

        //- Refine mesh
        mesh.update();

        #include "readTimeControls.H"

        // Warn if using too high of a courant number
        static bool hasWarned = false;
        if (maxCo > integrator.maxCo() && !hasWarned)
        {
            WarningInFunction
                << integrator.type() << " has a maximum stable Courant number "
                << "of " << integrator.maxCo() << " but a maximum Courant "
                << "number of " << maxCo << " has been specified" << endl;
            hasWarned = true;
        }

        #include "EigenCourantNos.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.name() << nl << endl;

        //- Move mesh
        mesh.move();

        //- Integrate the hyperbolic fluxes
        integrator.integrate();

        fluid.printInfo();

        integrator.clear();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
