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
    argList::addOption("elevation", "Elevation of the lowest point [m]");
    argList::addOption("phase", "Name of the phase");

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOobject atmospherePropertiesIO
    (
        "atmosphereProperties",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ
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
    wordList phases(1, phaseName);
    autoPtr<atmosphereModel> atmosphere;
    if (!atmospherePropertiesIO.typeHeaderOk<IOdictionary>(true))
    {
        scalar pRef = args.optionLookupOrDefault<scalar>
        (
            "pRef",
            Foam::constant::thermodynamic::Pstd
        );
        scalar elevation = args.optionLookupOrDefault<scalar>
        (
            "elevation",
            0
        );

        if (args.optionFound("phase"))
        {
            phases = args.optionRead<word>("phase");
        }
        else if (phaseProperties.found("phases"))
        {
            phases = phaseProperties.lookup<wordList>("phases");
        }
        atmosphere.set
        (
            new atmosphereModels::hydrostatic
            (
                mesh,
                pRef,
                elevation
            )
        );
    }
    else
    {
        IOdictionary atmosphereProperties(atmospherePropertiesIO);
        atmosphere =
            atmosphereModel::New(mesh, atmosphereProperties);

        if (atmosphereProperties.found("phase"))
        {
            phases = atmosphereProperties.lookup<word>("phase");
        }
        else if (phaseProperties.found("phases"))
        {
            phases = phaseProperties.lookup<wordList>("phases");
        }

    }

    autoPtr<fluidBlastThermo> thermo
    (
        fluidBlastThermo::New
        (
            phases.size(),
            mesh,
            phaseProperties
        )
    );

    Info<< "Initializing atmosphere." << endl;
    atmosphere->createAtmosphere(thermo());

    Info<< "Done" << nl << endl;
}


// ************************************************************************* //
