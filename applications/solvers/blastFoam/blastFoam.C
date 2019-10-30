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

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "staticFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "phaseCompressibleSystem.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;



    while (runTime.run())
    {
        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "estimateError.H"
        mesh.update();

        fluid->encode();

        Info<< "Calculating Fluxes" << endl;
        integrator->integrate();

        impulse +=
            (
                p
              - dimensionedScalar("pRef", dimPressure, 101298.0)
            )*runTime.deltaT();
        maxImpulse = max(maxImpulse, impulse);
        maxP = max(maxP, p);

        Info<< "max(p): " << max(p).value()
            << ", min(p): " << min(p).value() << endl;

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


        if(stopActive)
        {
            /* code */
            // check points in input file for vel mag and write/exit if > 0
            for(auto pI = 0; pI < stopPoints.size(); ++pI)
            {
                 auto UMagPoint = mag(U[mesh.findNearestCell(stopPoints[pI])]);
                 if(UMagPoint > stopEpsilon)
                 {
                    Info << "Stopping criteria detected at point: " << stopPoints[pI] << endl;
                    Info << "Mag velocity is : " << UMagPoint << endl;
                    Info<< "End\n" << endl;
                    runTime.writeAndEnd();
                    return 0;
                 }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
