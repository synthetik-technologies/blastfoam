/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

Application
    blastFoam

Description
    Multiphase compressible solver that uses Riemann solver to construct
    hyperbolic fluxes. Equation of states use the Mie–Grüneisen form.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "timeSelector.H"
#include "coupledMultiphaseCompressibleSystem.H"
#include "fvTimeIntegrator.H"

#include "parcelClouds.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    maxCo = min(maxCo, integrator.maxCo());
    scalar CoNum = fluid.CoNum();

    #include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        clouds.preUpdateMesh();
        integrator.preUpdateMesh();

        //- Refine the mesh
        mesh.update();

        scalar CoNum = fluid.CoNum();
        #include "readTimeControls.H"
        maxCo = min(maxCo, integrator.maxCo());

        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.name() << nl << endl;

        //- Move the mesh
        mesh.move();

        fluid.decode();
        clouds.evolve();
        alpha = clouds.alpha();

        fluid.eSource() = clouds.Sh(fluid.he());
        fluid.dragSource() = clouds.SU(fluid.U());

        Info<< "Calculating Fluxes" << endl;
        integrator.integrate();

        Info<< "max(p): " << max(p).value()
            << ", min(p): " << min(p).value() << endl;
        Info<< "max(T): " << max(T).value()
            << ", min(T): " << min(T).value() << endl;

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        integrator.clear();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
