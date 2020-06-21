/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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
    blastFoam

Description
    Multiphase compressible solver that uses Riemann solver to construct
    hyperbolic fluxes. Equation of states use the Mie–Grüneisen form.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "wedgeFvPatch.H"
#include "phaseCompressibleSystem.H"
#include "blastCompressibleTurbulenceModel.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "listTurbulenceModels",
        "List turbulenceModels"
    );
    #include "setRootCase.H"
    if (args.optionFound("listTurbulenceModels"))
    {
        Info<< "Turbulence models"
            << blast::laminarModel::dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "RAS models"
            << blast::RASModel::dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "LES models"
            << blast::LESModel::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        exit(1);
    }

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;



    while (runTime.run())
    {
        #include "eigenvalueCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

        fluid->encode();

        Info<< "Calculating Fluxes" << endl;
        integrator->integrate();

        Info<< "max(p): " << max(p).value()
            << ", min(p): " << min(p).value() << endl;
        Info<< "max(T): " << max(T).value()
            << ", min(T): " << min(T).value() << endl;

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        fluid->clearODEFields();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
