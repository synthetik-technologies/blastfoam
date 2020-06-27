/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate dynamic pressure
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

#include "dynamicPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dynamicPressure, 0);
    addToRunTimeSelectionTable(functionObject, dynamicPressure, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dynamicPressure::dynamicPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    rhoName_(dict.lookupOrDefault("rhoName", word("rho"))),
    UName_(dict.lookupOrDefault("UName", word("U"))),
    resultName_(dict.lookupOrDefault("name", word("dynamicP"))),
    store_(dict.lookupOrDefault("store", false))
{
    if (store_)
    {
        obr_.store
        (
            new volVectorField
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
                dimensionedVector("0", dimPressure, Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::dynamicPressure::~dynamicPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dynamicPressure::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("rhoName", rhoName_);
    dict.readIfPresent("UName", UName_);

    word origName = resultName_;
    bool origStore = store_;
    dict.readIfPresent("name", resultName_);
    dict.readIfPresent("store", store_);

    bool change = false;
    if ((origName != resultName_ && origStore) || (origStore && !store_))
    {
        change = true;
        clearObject(origName);
    }

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


bool Foam::functionObjects::dynamicPressure::execute()
{
    if
    (
        foundObject<volScalarField>(rhoName_)
     && foundObject<volVectorField>(UName_)
    )
    {
        const volScalarField& rho = lookupObject<volScalarField>(rhoName_);
        const volVectorField& U = lookupObject<volVectorField>(UName_);

        if (store_)
        {
            lookupObjectRef<volVectorField>(resultName_) =
                0.5*rho*mag(U)*U;
            return true;
        }

        return store
        (
            resultName_,
            0.5*rho*mag(U)*U
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


bool Foam::functionObjects::dynamicPressure::write()
{
    writeObject(resultName_);
    return true;
}


bool Foam::functionObjects::dynamicPressure::clear()
{
    if (!store_)
    {
        return clearObject(resultName_);
    }
    return true;
}

// ************************************************************************* //
