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
#include "thermophysicalTransportModel.H"
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
    intGradTName_
    (
        dict.lookupOrDefault
        (
            "intGradTName",
            IOobject::groupName("intGradT", IOobject::group(TName_))
        )
    ),
    intQExtName_
    (
        dict.lookupOrDefault
        (
            "intQExtName",
            IOobject::groupName("intQExt", IOobject::group(TName_))
        )
    ),
    intQName_
    (
        dict.lookupOrDefault
        (
            "intQName",
            IOobject::groupName("intQ", IOobject::group(TName_))
        )
    )
{
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

    const word thermophysicalTransportModelName
    (
        IOobject::groupName
        (
            thermophysicalTransportModel::typeName,
            IOobject::group(TName_)
        )
    );
    if
    (
        foundObject<thermophysicalTransportModel>
        (
            thermophysicalTransportModelName
        )
    )
    {
        if (!intQ_.valid())
        {
            intQ_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        intQName_,
                        Time::timeName(time_.startTime().value()),
                        mesh_,
                        restartOnRestart_
                      ? IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimMass/sqr(dimTime), 0.0)
                )
            );
        }
        const thermophysicalTransportModel& ttm =
            lookupObject<thermophysicalTransportModel>
            (
                thermophysicalTransportModelName
            );
        intQ_() = intQ_->oldTime() + deltaT*fvc::average(ttm.q());
    }
    else
    {
        if (!intGradT_.valid())
        {
            intGradT_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        intGradTName_,
                        Time::timeName(time_.startTime().value()),
                        mesh_,
                        restartOnRestart_
                      ? IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimTemperature/dimLength*dimTime, 0.0)
                )
            );
        }

        intGradT_() =
            intGradT_->oldTime()
          + fvc::average
            (
                fvc::snGrad(mesh_.lookupObject<volScalarField>(TName_))
            )*deltaT;
    }

    if (mesh_.foundObject<volScalarField>(qrName_))
    {
        if (!intQExt_.valid())
        {
            intQExt_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        intQExtName_,
                        Time::timeName(time_.startTime().value()),
                        mesh_,
                        restartOnRestart_
                      ? IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimMass/sqr(dimTime), 0.0)
                )
            );
        }

        const volScalarField& qr =
            mesh_.lookupObject<volScalarField>(qrName_);
        intQExt_() = intQExt_->oldTime() + deltaT*qr;

        if (intQ_.valid())
        {
            intQ_() += deltaT*qr;
        }
    }



    return true;
}


bool Foam::functionObjects::thermalImpulse::write()
{
    if (obr_.time().timeIndex() == obr_.time().startTimeIndex())
    {
        return true;
    }
    bool good = true;

    if (intGradT_.valid())
    {
        good = intGradT_->write() && good;
    }
    if (intQExt_.valid())
    {
        good = intQExt_->write() && good;
    }
    if (intQ_.valid())
    {
        good = intQ_->write() && good;
    }
    return good;
}


// ************************************************************************* //
