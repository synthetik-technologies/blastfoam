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
#include "fvc.H"
#include "fluidBlastThermo.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

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
            IOobject::MUST_READ
        )
    );

    autoPtr<atmosphereModel> atmosphere
    (
        atmosphereModel::New(mesh, atmosphereProperties)
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

    //- Name of phase (Default = word::null)
    word phaseName
    (
        atmosphereProperties.lookupOrDefault("phaseName", word::null)
    );

    //- Shared switches
    Switch sharedPressure
    (
        atmosphereProperties.lookupOrDefault("sharedPressure", true)
    );
    Switch sharedTemperature
    (
        atmosphereProperties.lookupOrDefault("sharedTemperature", true)
    );

    wordList phases(phaseProperties.lookupOrDefault("phases", wordList(1,word::null)));
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

    forAll(thermo->p().boundaryField(), patchi)
    {
        thermo->p().boundaryFieldRef()[patchi] =
            thermo->p().boundaryField()[patchi].patchInternalField();
    }
    thermo->T().write();
    thermo->rho().write();
    thermo->p().write();

    Info<< "Done" << nl << endl;
}


// ************************************************************************* //
