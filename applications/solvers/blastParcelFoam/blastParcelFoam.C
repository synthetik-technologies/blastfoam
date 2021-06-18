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
#include "zeroGradientFvPatchFields.H"
#include "wedgeFvPatch.H"
#include "coupledMultiphaseCompressibleSystem.H"
#include "blastFluidThermoMomentumTransportModel.H"
#include "timeIntegrator.H"
#include "fluidThermoModel.H"

#include "basicThermoCloud.H"

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
        parcels.storeGlobalPositions();

        //- Refine the mesh
        mesh.refine();

        #include "eigenvalueCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //- Move the mesh
        mesh.update();

        fluid.decode();
        parcels.evolve();
        theta = parcels.theta();
        fluid.encode();

        fluid.eSource() = parcels.Sh(fluid.e());
        fluid.dragSource() = parcels.SU(fluid.U());

        Info<< "Calculating Fluxes" << endl;
        integrator->integrate();

        //- Clear the flux scheme
        fluid.flux().clear();

        Info<< "max(p): " << max(p).value()
            << ", min(p): " << min(p).value() << endl;
        Info<< "max(T): " << max(T).value()
            << ", min(T): " << min(T).value() << endl;

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        integrator->clearODEFields();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
