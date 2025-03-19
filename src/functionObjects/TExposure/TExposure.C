/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate TExposure
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

#include "TExposure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TExposure, 0);
    addToRunTimeSelectionTable(functionObject, TExposure, dictionary);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::TExposure::burnModel,
    2
>::names[] = {"CEM43", "Henriques"};

const Foam::NamedEnum
<
    Foam::functionObjects::TExposure::burnModel,
    2
> Foam::functionObjects::TExposure::burnModelNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TExposure::TExposure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    TName_(dict.lookupOrDefault("TName", word("T"))),
    exposure_
    (
        IOobject
        (
            dict.lookupOrDefault
            (
                "fieldName",
                IOobject::groupName("TExposure", IOobject::group(TName_))
            ),
            runTime.timeName(),
            mesh_,
            restartOnRestart_
          ? IOobject::NO_READ
          : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimTime, 0.0)
    ),
    burnModel_(burnModel::CEM43),
    P_(0.0, 0.0),
    EbyR_(1.0, 1.0),
    Tswitch_(0.0),
    Tmin_(0.0),
    Tmax_(273.15 + 300.0)
{
    executeAtStart_ = false;
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::TExposure::~TExposure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TExposure::read(const dictionary& dict)
{
    Log << type() << " " << name() << ":" << nl;
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("restartOnRestart", restartOnRestart_);

    burnModel_ = burnModelNames_[dict.lookup<word>("burnModel")];
    if (burnModel_ == burnModel::HENRIQUES)
    {
        word source("Henriques");
        dict.readIfPresent("source", source);
        if (source == "Henriques")
        {
            P_ = {3.1e98, 3.1e98};
            EbyR_ = {75000, 75000};
        }
        else if (source == "WeaverStoll")
        {
            P_ = {2.185e124, 1.823e51};
            EbyR_ = {93534.9, 39109.8};
            Tswitch_ = 50.0 + 273.15;
        }
        else if (source == "Takata")
        {
            P_ = {4.32e64, 9.39e104};
            EbyR_ = {50000.0, 80000.0};
            Tswitch_ = 50.0 + 273.15;
            Tmin_ = 44.0 + 273.15;
            Tmax_ = 60.0 + 273.15;
        }
        else if (source == "MehtaWong")
        {
            EbyR_ = {55000, 55000};
            word skinLayer(dict.lookup("skinLayer"));
            if (skinLayer == "epidermis")
            {
                P_ = {1.43e72, 1.43e72};
            }
            else if (skinLayer == "dermis")
            {
                P_ = {2.86e69, 2.86e69};
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Unknown skinLayer " << skinLayer << nl
                    << "Valid skinLayers are:" << nl
                    << "    epidermis" << nl
                    << "    dermis" << nl
                    << endl
                    << abort(FatalIOError);
            }
        }
        else if (source == "custom")
        {
            dict.lookup("preExponential") >> P_;
            dict.lookup("deltaEbyR") >> EbyR_;
            if (P_.first() != P_.second() && EbyR_.first() != EbyR_.second())
            {
                dict.lookup("TSwitch") >> Tswitch_;
            }
            dict.readIfPresent("Tmin", Tmin_);
            dict.readIfPresent("Tmax", Tmax_);
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Unknown epidermis source " << source
                << "Valid sources are: " << nl
                << "    Henriques" << nl
                << "    WeaverStoll" << nl
                << "    Takata" << nl
                << "    MehtaWong" << nl
                << "    custom" << nl
                << endl
                << abort(FatalIOError);
        }
    }
    else
    {
        Tswitch_ = 43.0 + 273.15;
    }


    Log << endl;

    return true;
}


bool Foam::functionObjects::TExposure::execute()
{
    const volScalarField& T =
        mesh_.lookupObject<volScalarField>(TName_);

    const scalar dt = obr_.time().deltaTValue();
    const volScalarField& exposure0 = exposure_.oldTime();
    volScalarField::Boundary& bexposure = exposure_.boundaryFieldRef();

    bool warn = false;
    if (burnModel_ == burnModel::CEM43)
    {
        forAll(T, celli)
        {
            scalar Ti = T[celli];
            if (Ti > Tmin_)
            {
                if (Ti > Tmax_)
                {
                    warn = true;
                    Ti = Tmax_;
                }
                exposure_[celli] =
                    exposure0[celli]
                  + (
                        Ti > Tswitch_
                      ? pow(0.5, Tswitch_ - Ti)
                      : pow(0.25, Tswitch_ - Ti)
                    )*dt;
            }
        }

        forAll(bexposure, patchi)
        {
            const scalarField& pT = T.boundaryField()[patchi];
            const scalarField& pexposure0 = exposure0.boundaryField()[patchi];
            scalarField& pexposure = bexposure[patchi];
            forAll(pexposure, facei)
            {
                scalar Ti = pT[facei];
                if (Ti > Tmin_)
                {
                    if (Ti > Tmax_)
                    {
                        warn = true;
                        Ti = Tmax_;
                    }
                    pexposure[facei] =
                        pexposure0[facei]
                      + (
                            Ti > Tswitch_
                          ? pow(0.5, Tswitch_ - Ti)
                          : pow(0.25, Tswitch_ - Ti)
                        )*dt;
                }
            }
        }
    }
    else
    {
        forAll(T, celli)
        {
            scalar Ti = T[celli];
            if (Ti > Tmin_)
            {
                if (Ti > Tmax_)
                {
                    warn = true;
                    Ti = Tmax_;
                }
                exposure_[celli] =
                    exposure0[celli]
                  + (
                        Ti > Tswitch_
                      ? P_.second()*exp(-EbyR_.second()/Ti)
                      : P_.first()*exp(-EbyR_.first()/Ti)
                    )*dt;
            }
        }

        forAll(bexposure, patchi)
        {
            const scalarField& pT = T.boundaryField()[patchi];
            const scalarField& pexposure0 = exposure0.boundaryField()[patchi];
            scalarField& pexposure = bexposure[patchi];
            forAll(pexposure, facei)
            {
                scalar Ti = pT[facei];
                if (Ti > Tmin_)
                {
                    if (Ti > Tmax_)
                    {
                        warn = true;
                        Ti = Tmax_;
                    }
                    pexposure[facei] =
                        pexposure0[facei]
                      + (
                            Ti > Tswitch_
                          ? P_.second()*exp(-EbyR_.second()/Ti)
                          : P_.first()*exp(-EbyR_.first()/Ti)
                        )*dt;
                }

            }
        }
    }

    if (debug && warn)
    {
        WarningInFunction
            << "Temperature surpassed maximum model temperature, "
            << Tmax_ << endl;
    }

    return true;
}


bool Foam::functionObjects::TExposure::write()
{
    if (obr_.time().timeIndex() == obr_.time().startTimeIndex())
    {
        return true;
    }
    return exposure_.write();
}


// ************************************************************************* //
