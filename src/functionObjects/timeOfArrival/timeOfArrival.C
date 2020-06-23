/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
20-06-2020 Jeff Heylmun:    | Time of arrival implementation
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

#include "timeOfArrival.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeOfArrival, 0);
    addToRunTimeSelectionTable(functionObject, timeOfArrival, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::volScalarField&
Foam::functionObjects::timeOfArrival::lookupOrCreate
(
    const word& name,
    const dimensionSet& dims
) const
{
    Log << "    Reading/initialising field " << name << endl;

    if (obr_.foundObject<volScalarField>(name))
    {
        return obr_.lookupObjectRef<volScalarField>(name);
    }

    // Store on registry
    volScalarField* fieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                name,
                obr_.time().timeName(),
                obr_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("0", dims, 0.0)
        )
    );
    fieldPtr->store(fieldPtr);

    return *fieldPtr;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeOfArrival::timeOfArrival
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pMax_
    (
        lookupOrCreate
        (
            IOobject::groupName
            (
                IOobject::member(pName_)
              + "Max",
                IOobject::group(pName_)
            ),
            dimPressure
        )
    ),
    timeOfArrival_
    (
        lookupOrCreate
        (
            IOobject::groupName
            (
                "timeOfArrival",
                IOobject::group(pName_)
            ),
            dimTime
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::timeOfArrival::~timeOfArrival()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeOfArrival::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::timeOfArrival::execute()
{
    const fvMesh& mesh = this->mesh_;
    const volScalarField& p(mesh.lookupObject<volScalarField>(pName_));

    forAll(pMax_, celli)
    {
        if (p[celli] > pMax_[celli])
        {
            timeOfArrival_[celli] = obr_.time().value();
            pMax_[celli] = p[celli];
        }
    }
    forAll(pMax_.boundaryField(), patchi)
    {
        forAll(pMax_.boundaryField()[patchi], facei)
        {
            if
            (
                p.boundaryField()[patchi][facei]
              > pMax_.boundaryField()[patchi][facei]
            )
            {
                timeOfArrival_.boundaryFieldRef()[patchi][facei] =
                    obr_.time().value();
                pMax_.boundaryFieldRef()[patchi][facei] =
                    p.boundaryField()[patchi][facei];
            }
        }
    }

    return true;
}


bool Foam::functionObjects::timeOfArrival::write()
{
    pMax_.write();
    timeOfArrival_.write();
    return true;
}


// ************************************************************************* //
