/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate thermalImpulse
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "thermalImpulse.H"
#include "fvcAverage.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(thermalImpulse, 0);
    addToRunTimeSelectionTable(functionObject, thermalImpulse, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::thermalImpulse::thermalImpulse
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    TName_(dict.lookupOrDefault("TName", word("T"))),
    qrName_(dict.lookupOrDefault("qrName", word("qr"))),
    intGradT_
    (
        IOobject
        (
            dict.lookupOrDefault
            (
                "intGradTName",
                IOobject::groupName("intGradT", IOobject::group(TName_))
            ),
            runTime.timeName(),
            mesh_,
            restartOnRestart_
          ? IOobject::NO_READ
          : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimTemperature/dimLength*dimTime, 0.0)
    ),
    intQExt_
    (
        IOobject
        (
            dict.lookupOrDefault
            (
                "intQName",
                IOobject::groupName("intQExt", IOobject::group(TName_))
            ),
            runTime.timeName(),
            mesh_,
            restartOnRestart_
          ? IOobject::NO_READ
          : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimMass/sqr(dimTime), 0.0)
    )
{
    executeAtStart_ = false;
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::thermalImpulse::~thermalImpulse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::thermalImpulse::read(const dictionary& dict)
{
    Log << type() << " " << name() << ":" << nl;
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("restartOnRestart", restartOnRestart_);

    Log << endl;

    return true;
}


bool Foam::functionObjects::thermalImpulse::execute()
{
    const dimensionedScalar& deltaT = mesh_.time().deltaT();

    intGradT_ =
        intGradT_.oldTime()
      + deltaT
       *fvc::average
        (
            fvc::snGrad(mesh_.lookupObject<volScalarField>(TName_))
        );
    if (mesh_.foundObject<volScalarField>(qrName_))
    {
        intQExt_ =
            intQExt_.oldTime()
          + deltaT*mesh_.lookupObject<volScalarField>(qrName_);
    }

    return true;
}


bool Foam::functionObjects::thermalImpulse::write()
{
    return intGradT_.write() && intQExt_.write();
}


// ************************************************************************* //
