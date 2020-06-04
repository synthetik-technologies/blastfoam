/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
29-4-2019 Jeff Heylmun:    Added moment transort
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
    oneWayCoupledHyperbolicFoam

Description
    Moment transport of velocity based NDF that interacts with a constant
    velocity field by means of drag.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "staticFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "phaseCompressibleSystem.H"
#include "timeIntegrator.H"
#include "fluidThermoModel.H"
#include "populationBalanceModel.H"
#include "quadratureApproximations.H"
#include "mappedPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    #include "CourantNos.H"
    #include "setInitialDeltaT.H"



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "CourantNos.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

        fluid->encode();

        Info<< "Calculating Fluxes" << endl;
        integrator->integrate();

        Info<< "Solving PBE" << endl;
        populationBalance->solve();

        Info<< "Solving velocity abscissae" << endl;
        #include "vEqns.H"

        Info<< "Calculating mean quantities"<< endl;
        #include "computeParticleFields.H"

        fluid->clearODEFields();

        Info<< "max(p): " << max(p).value()
            << ", min(p): " << min(p).value() << endl;

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
