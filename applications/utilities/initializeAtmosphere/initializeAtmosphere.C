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
    Initialize an atmosphere using either the U.S. standard 76 model, standard
    atmosphere from NASA, or using a simple hydrostatic pressure.

    References:
    \verbatim
        The U.S Standard Atmosphere, 1976. U.S. Government Printing Office.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "HashSet.H"
#include "fvc.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fluidThermoModel.H"
#include "lookupTable1D.H"
#include "twoPhaseFluidThermo.H"
#include "multiphaseFluidThermo.H"
#include "PtrListDictionary.H"

#include "atmosphereFuncs.H"

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
            runTime.system(),
            mesh,
            IOobject::MUST_READ
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

    //- Name of phase (Default = word::null)
    word phaseName(atmosphereProperties.lookupOrDefault("phaseName", word::null));

    //- Shared switches
    Switch sharedPressure(atmosphereProperties.lookupOrDefault("sharedPressure", true));
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

    HashSet<> atmosphereTypes;
    atmosphereTypes.insert("simple");
    atmosphereTypes.insert("tabulated");
    atmosphereTypes.insert("standard76");

    word model(atmosphereProperties.lookup("type"));

    dimensionedScalar groundElevation
    (
        "groundElevation",
        dimLength,
        atmosphereProperties
    );

    Info<< "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
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
    fluidThermoModel* thermoPtr;
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
    fluidThermoModel& thermo(*thermoPtr);

    vector dir((-g/mag(g)).value());
    volScalarField h("h", dir & mesh.C());
    h += (groundElevation - min(h));

    Info<< "Initializing " << model <<" atmosphere." << endl;
    if (model == "simple")
    {
        //- Simple hydrostatic equilibrium
        simpleAtmosphere(atmosphereProperties, g, dir, h, p, rho);
    }
    else if (model == "tabulated")
    {
        //- Atmosphere is based on a lookup table of pressure and temperature Vs.
        //  height. Air is the hard coded
        tableAtmosphere(atmosphereProperties, g, dir, h, p, rho);
    }
    else if (model == "standard76")
    {
        //- U.S. standard 76 atmosphere model
        standardAtmosphere76(atmosphereProperties, g, dir, h, p, rho);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown " << model << " atmosphere model " << nl
            << "Valid atmosphere model types are:"
            << atmosphereTypes << nl
            <<abort(FatalError);
    };
    e = thermo.calce();
    T = thermo.calcT();

    if (atmosphereProperties.found("equilibriumTemperature") && multiphase)
    {
        forAll(phases, phasei)
        {
            thermos[phasei].T() = T;
            thermos[phasei].T().write();
        }
    }

    p.write();
    rho.write();
    T.write();

    Info<< "Done" << nl << endl;
}


// ************************************************************************* //
