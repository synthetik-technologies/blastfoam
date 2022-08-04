/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate impulse
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

#include "impulse.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(impulse, 0);
    addToRunTimeSelectionTable(functionObject, impulse, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::impulse::impulse
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_("pRef", dimPressure, dict),
    impulse_
    (
        IOobject
        (
            IOobject::groupName("impulse", IOobject::group(pName_)),
            runTime.timeName(),
            mesh_,
            restartOnRestart_
          ? IOobject::NO_READ
          : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimPressure*dimTime, 0.0)
    )
{
    if (!dict.lookupOrDefault("executeAtStart", false))
    {
        executeAtStart_ = false;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::impulse::~impulse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::impulse::read(const dictionary& dict)
{
    Log << type() << " " << name() << ":" << nl;
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    pRef_.read(dict);

    Log << endl;

    return true;
}


bool Foam::functionObjects::impulse::execute()
{
    const volScalarField& p =
        mesh_.lookupObject<volScalarField>(pName_);

    impulse_ = impulse_.oldTime() + (p - pRef_)*obr_.time().deltaT();

    return true;
}


bool Foam::functionObjects::impulse::write()
{
    return impulse_.write();
}


// ************************************************************************* //
