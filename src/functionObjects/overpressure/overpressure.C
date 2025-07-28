/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies: | Calculate overpressure
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
    read(dict);
    if (store_)
    {
        obr_.store
        (
            new volScalarField
            (
                IOobject
                (
                    resultName_,
                    obr_.time().name(),
                    obr_
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
    if (!p0Ptr_.valid() && dict.lookupOrDefault("nonUniformPRef", false))
    {
        typeIOobject<volScalarField> p0IO
        (
            IOobject::groupName("p0", IOobject::group(pName_)),
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if (p0IO.headerOk())
        {
            p0Ptr_.set(new volScalarField(p0IO, mesh_));
        }
        else
        {
            p0IO.readOpt() = IOobject::NO_READ;
            p0Ptr_.set
            (
                new volScalarField
                (
                    p0IO,
                    mesh_.lookupObject<volScalarField>(pName_)
                )
            );
        }
    }

    if (!p0Ptr_.valid())
    {
        pRef_.read(dict);
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
            if (p0Ptr_.valid())
            {
                lookupObjectRef<volScalarField>(resultName_) = p - p0Ptr_();
            }
            else
            {
                lookupObjectRef<volScalarField>(resultName_) = p - pRef_;
            }
            return true;
        }

        if (p0Ptr_.valid())
        {
            return store(resultName_, p - p0Ptr_());
        }
        else
        {
            return store(resultName_, p - pRef_);
        }
    }
    else
    {
        WarningInFunction
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        return false;
    }
}


bool Foam::functionObjects::overpressure::write()
{
    return writeObject(resultName_);
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
