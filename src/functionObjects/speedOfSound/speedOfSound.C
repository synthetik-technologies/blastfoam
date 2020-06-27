/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate speed of sound with blastFoam thermo
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

#include "speedOfSound.H"
#include "fluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(speedOfSound, 0);
    addToRunTimeSelectionTable(functionObject, speedOfSound, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::speedOfSound::speedOfSound
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookupOrDefault("phaseName", word::null)),
    resultName_(IOobject::groupName("speedOfSound", phaseName_)),
    store_(dict.lookupOrDefault("store", false))
{
    if (store_)
    {
        obr_.store
        (
            new volScalarField
            (
                IOobject
                (
                    resultName_,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("0", dimVelocity, Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::speedOfSound::~speedOfSound()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::speedOfSound::read
(
    const dictionary& dict
)
{
    word origName = phaseName_;
    bool origStore = store_;
    dict.readIfPresent("phaseName", phaseName_);
    dict.readIfPresent("store", store_);

    bool change = false;
    if ((origName != phaseName_ && origStore) || (origStore && !store_))
    {
        change = true;
        clearObject(resultName_);
    }

    resultName_ = IOobject::groupName("speedOfSound", phaseName_);
    if (store_ && change)
    {
        obr_.store
        (
            new volScalarField
            (
                IOobject
                (
                    resultName_,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("0", dimVelocity, Zero)
            )
        );
    }

    return true;
}


bool Foam::functionObjects::speedOfSound::execute()
{
    if
    (
        foundObject<fluidThermoModel>
        (
            IOobject::groupName("basicThermo", phaseName_)
        )
    )
    {
        tmp<volScalarField> c =
            lookupObject<fluidThermoModel>
            (
                IOobject::groupName
                (
                    "basicThermo",
                    phaseName_
                )
            ).speedOfSound();
        if (store_)
        {
            lookupObjectRef<volScalarField>(resultName_) = c;
            return true;
        }

        return store
        (
            resultName_,
            c
        );
    }
    else
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        return false;
    }
}


bool Foam::functionObjects::speedOfSound::write()
{
    writeObject(resultName_);
    return true;
}


bool Foam::functionObjects::speedOfSound::clear()
{
    if (!store_)
    {
        return clearObject(resultName_);
    }
    return true;
}

// ************************************************************************* //
