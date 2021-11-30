/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    elasticSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic solid bodies.

    Displacement field U is solved for using a total Lagrangian approach,
    also generating the strain tensor field epsilon and stress tensor
    field sigma.

    With optional multi-material solid interface correction ensuring
    correct tractions on multi-material interfaces

Author
    Philip Cardiff
    multi-material by Tukovic et al. 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicBlastFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    {
        scalar dTOld(runTime.deltaTValue());
        solid.setDeltaT(runTime);
        runTime.setDeltaT(min(runTime.deltaTValue(), dTOld));
    }
    scalar CoNum = solid.CoNum();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        CoNum = solid.CoNum();
        maxCo = min(maxCo, solid.maxCoNum());
        Info<< "Max Courant Number = " << CoNum << endl;
        #include "setDeltaT.H"
//         solid.setDeltaT(runTime);
//         #include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        solid.evolve();

        solid.updateTotalFields();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl << endl;

    }

    solid.end();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
