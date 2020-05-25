/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate overpressure
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

#include "overpressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(overpressure, 0);
    addToRunTimeSelectionTable(functionObject, overpressure, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::overpressure::overpressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    resultName_(IOobject::groupName("overpressure", IOobject::group(pName_))),
    pRef_("pRef", dimPressure, dict),
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
                dimensionedScalar("0", dimPressure, Zero)
            )
        );
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::overpressure::~overpressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::overpressure::read
(
    const dictionary& dict
)
{
    word origName = pName_;
    bool origStore = store_;
    dict.readIfPresent("pName", pName_);
    dict.readIfPresent("store", store_);

    bool change = false;
    if ((origName != pName_ && origStore) || (origStore && !store_))
    {
        change = true;
        clearObject(resultName_);
    }

    resultName_ = IOobject::groupName("overpressure", IOobject::group(pName_));
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
                dimensionedScalar("0", dimPressure, Zero)
            )
        );
    }

    return true;
}


bool Foam::functionObjects::overpressure::execute()
{
    if (foundObject<volScalarField>(pName_))
    {
        const volScalarField& p(lookupObject<volScalarField>(pName_));

        if (store_)
        {
            lookupObjectRef<volScalarField>(resultName_) = p - pRef_;
            return true;
        }

        return store
        (
            resultName_,
            p - pRef_
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


bool Foam::functionObjects::overpressure::write()
{
    writeObject(resultName_);
    return true;
}


bool Foam::functionObjects::overpressure::clear()
{
    if (!store_)
    {
        return clearObject(resultName_);
    }
    return true;
}

// ************************************************************************* //
