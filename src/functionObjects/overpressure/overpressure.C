/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    resultName_(IOobject::groupName("overPressure", IOobject::group(pName_))),
    pRef_("pRef", dimPressure, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::overpressure::~overpressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::overpressure::read
(
    const dictionary& dict
)
{
    pRef_.read(dict);
    pName_ = dict.lookupOrDefault("pName", word("p"));
    return true;
}


bool Foam::functionObjects::overpressure::execute()
{
    if (foundObject<volScalarField>(pName_))
    {
        const volScalarField& p(lookupObject<volScalarField>(pName_));

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
    return clearObject(resultName_);
}

// ************************************************************************* //
