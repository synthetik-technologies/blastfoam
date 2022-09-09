/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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

Description
    Initializes a pressure and density field for gravitational stability.

    References:
    \verbatim
        The U.S Standard Atmosphere, 1976. U.S. Government Printing Office.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "atmosphereModel.H"
#include "hydrostaticAtmosphereModel.H"
#include "fvc.H"
#include "fluidBlastThermo.H"
#include "thermodynamicConstants.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("pRef", "Reference pressure [Pa]");
    argList::addOption("hRef", "Height of the lowest point [m]");
    argList::addOption("phase", "Name of the phase");
    argList::addOption("zone", "Cell zone to set");
    argList::addOption("type", "Model used to set the fields");
    argList::addOption("fixedPatches", "patches to fix pressure on");
    argList::addOption("refCell", "Reference cell");
    argList::addOption("refPoint", "Reference point");

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary atmosphereProperties
    (
        IOobject
        (
            "atmosphereProperties",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT
        )
    );
    word type = atmosphereModels::hydrostatic::typeName;
    if (atmosphereProperties.headerOk())
    {
        type = atmosphereProperties.lookup<word>("type");
    }
    else if (args.optionFound("type"))
    {
        type = args.optionRead<word>("type");
    }

    if (args.optionFound("hRef"))
    {
        atmosphereProperties.set("hRef", args.optionRead<scalar>("hRef"));
    }
    else if (!atmosphereProperties.found("hRef"))
    {
        WarningInFunction << "hRef was not provided, using 0" << endl;
    }

    label refSet = 0;
    if (args.optionFound("fixedPatches"))
    {
        refSet = 2;
        atmosphereProperties.set
        (
            "fixedPatches",
            args.optionRead<wordReList>("fixedPatches")
        );
    }

    label zoneID = -1;
    if (args.optionFound("zone"))
    {
        zoneID = mesh.cellZones()[args.optionRead<word>("zone")].index();
    }
    autoPtr<atmosphereModel> atmosphere
    (
        atmosphereModel::New
        (
            type,
            mesh,
            atmosphereProperties,
            zoneID
        )
    );

    IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    word phaseName = word::null;
    wordList phases(1, word::null);
    if (args.optionFound("phase"))
    {
        phases = args.optionRead<word>("phase");
        if (phaseProperties.found("phases"))
        {
            phaseName = phases[0];
        }
    }
    else if (atmosphereProperties.found("phase"))
    {
        phases = atmosphereProperties.lookup<word>("phase");
    }
    else if (phaseProperties.found("phases"))
    {
        phases = phaseProperties.lookup<wordList>("phases");
    }

    autoPtr<fluidBlastThermo> thermo
    (
        fluidBlastThermo::New
        (
            phases.size(),
            mesh,
            phaseProperties,
            phaseName
        )
    );

    if (thermo->p().needReference())
    {
        if (args.optionFound("pRef"))
        {
            atmosphereProperties.set("pRefValue", args.optionRead<scalar>("pRef"));
            refSet++;
        }
        else if (atmosphereProperties.found("pRef"))
        {
            refSet++;
        }

        if (args.optionFound("refCell"))
        {
            atmosphereProperties.set("pRefCell", args.optionRead<int>("refCell"));
            refSet++;
        }
        else if (args.optionFound("refPoint"))
        {
            atmosphereProperties.set("pRefPoint", args.optionRead<vector>("refPoint"));
            refSet++;
        }
        else if
        (
            atmosphereProperties.found("pRef")
         || atmosphereProperties.found("pRefCell")
        )
        {
            refSet++;
        }

        if (refSet < 2)
        {
            FatalErrorInFunction
                << "Could not determine reference pressure state" << nl
                << "please provide pRef and refCell/pRefCell or refPoint/pRefPoint" << nl
                << " or provide fixed pressure patches" << endl
                << abort(FatalError);
        }
    }

    Info<< "Initializing atmosphere." << endl;
    atmosphere->createAtmosphere(thermo());

    runTime.writeNow();

    Info<< "Done" << nl << endl;
}


// ************************************************************************* //
