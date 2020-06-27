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
#include "fluidThermoModel.H"
#include "twoPhaseFluidThermo.H"
#include "multiphaseFluidThermo.H"
#include "PtrListDictionary.H"


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

    //- Field names
    word rhoName(IOobject::groupName("rho", phaseName));
    word pName
    (
        sharedPressure
      ? "p"
      : IOobject::groupName("p", phaseName)
    );
    word TName
    (
        sharedTemperature
      ? "T"
      : IOobject::groupName("T", phaseName)
    );
    word eName
    (
        sharedTemperature
      ? "e"
      : IOobject::groupName("e", phaseName)
    );

    //- Density field
    volScalarField rho
    (
        IOobject
        (
            rhoName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField p
    (
        IOobject
        (
            pName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField T
    (
        IOobject
        (
            TName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField e
    (
        IOobject
        (
            eName,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimVelocity), -1.0),
        fluidThermoModel::eBoundaryTypes(T),
        fluidThermoModel::eBoundaryBaseTypes(T)
    );

    wordList phases(phaseProperties.lookupOrDefault("phases", wordList()));
    bool multiphase(false);
    if (phases.size() > 0)
    {
        multiphase = phaseProperties.subDict(phases[0]).found("thermoCoeffs");
    }
    PtrListDictionary<fluidThermoModel> thermos
    (
        multiphase
      ? phases.size()
      : 1
    );
    fluidThermoModel* thermoPtr = nullptr;
    if (!multiphase)
    {
        if (phases.size() <= 1)
        {
            thermos.set
            (
                0,
                "mixture",
                fluidThermoModel::New
                (
                    word::null,
                    p,
                    rho,
                    e,
                    T,
                    phaseProperties.subDict("mixture"),
                    true
                ).ptr()
            );
            thermoPtr = &thermos[0];
        }
        else if (phases.size() == 2)
        {
            thermos.set
            (
                0,
                "mixture",
                twoPhaseFluidThermo::New
                (
                    word::null,
                    p,
                    rho,
                    e,
                    T,
                    phaseProperties,
                    true
                ).ptr()
            );
            thermoPtr =
                thermos[0].thermo(0).rho().group() == phaseName
              ? &thermos[0].thermo(0)
              : &thermos[0].thermo(1);
        }
        else
        {
            thermos.set
            (
                0,
                "mixture",
                multiphaseFluidThermo::New
                (
                    word::null,
                    p,
                    rho,
                    e,
                    T,
                    phaseProperties,
                    true
                ).ptr()
            );
            forAll(phases, phasei)
            {
                if (thermos[0].thermo(phasei).rho().group() == phaseName)
                {
                    thermoPtr = &thermos[0].thermo(phasei);
                    break;
                }
            }
        }
    }
    else
    {
        forAll(phases, phasei)
        {
            thermos.set
            (
                phasei,
                phases[phasei],
                fluidThermoModel::New
                (
                    phaseName,
                    p,
                    rho,
                    e,
                    T,
                    phaseProperties.subDict(phases[phasei]).subDict("thermoCoeffs"),
                    true
                ).ptr()
            );
            if (phases[phasei] == phaseName)
            {
                thermoPtr = &thermos[phasei];
            }
        }
    }

    Info<< "Initializing atmosphere." << endl;
    atmosphere->createAtmosphere(p, rho);

    if (thermoPtr)
    {
        e = thermoPtr->calce();
        T = thermoPtr->calcT();
    }

    if (multiphase)
    {
        if (atmosphereProperties.lookupType<Switch>("equilibriumTemperature"))
        {
            forAll(phases, phasei)
            {
                Info<< "Setting " << thermos[phasei].T().name() << endl;
                thermos[phasei].T() = T;
                thermos[phasei].T().write();
            }
        }
    }

    p.write();
    rho.write();
    T.write();

    Info<< "Done" << nl << endl;
}


// ************************************************************************* //
