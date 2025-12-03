/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
15-01-2020 Jeff Heylmun     : Added Riemann based fluxes
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    blastXiFoam

Description
    Solver for highly compressible premixed/partially-premixed
    combustion with turbulence modeling.

    Combusting RANS code using the b-Xi two-equation model.
    Xi may be obtained by either the solution of the Xi transport
    equation or from an algebraic expression.  Both approaches are
    based on Gulder's flame speed correlation which has been shown
    to be appropriate by comparison with the results from the
    spectral model.

    Strain effects are encorporated directly into the Xi equation
    but not in the algebraic approximation.  Further work need to be
    done on this issue, particularly regarding the enhanced removal rate
    caused by flame compression.  Analysis using results of the spectral
    model will be required.

    For cases involving very lean Propane flames or other flames which are
    very strain-sensitive, a transport equation for the laminar flame
    speed is present.  This equation is derived using heuristic arguments
    involving the strain time scale and the strain-rate at extinction.
    the transport velocity is the same as that for the Xi equation.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "psiuCompressibleSystem.H"
#include "fluidThermophysicalTransportModel.H"
#include "fluxScheme.H"
#include "fvTimeIntegrator.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    fluxSchemeBase::needEnergyFlux = true;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    maxCo = min(maxCo, integrator.maxCo());
    scalar CoNum = fluid.CoNum();
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // Update fvModels and constraints
        integrator.preUpdateMesh();

        // Refine the mesh
        mesh.update();

        // Update Courant number
        CoNum = fluid.CoNum();
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

        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.name() << nl << endl;

        Info<< "Calculating Fluxes" << endl;
        integrator.integrate();

//         #include "ftEqn.H"
//         #include "bEqn.H"

        Info<< "    max(p) = " << max(p).value()
            << ", min(p) = " << min(p).value() << nl
            << "    max(T) = " << max(T).value()
            << ", min(T) = " << min(T).value() << nl
            << "    max(b) = " << max(b).value()
            << ", min(b) = " << min(b).value() << nl
            << "    Combustion progress = "
            << 100*(scalar(1) - b)().weightedAverage(mesh.V()).value() << "%"
            << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
