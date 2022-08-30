/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies: | Calculate dynamic pressure
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::dynamicPressure::calc()
{
    bool good = true;
    if (!foundObject<volScalarField>(rhoName_))
    {
        cannotFindObject(rhoName_);
        good = false;
    }
    if (!foundObject<volVectorField>(UName_))
    {
        cannotFindObject(UName_);
        good = false;
    }

    if (!good)
    {
        return false;
    }

    const volScalarField& rho = lookupObject<volScalarField>(rhoName_);
    const volVectorField& U = lookupObject<volVectorField>(UName_);

    return store
    (
        resultName_,
        0.5*rho*mag(U)*U
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dynamicPressure::dynamicPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "dynamicP"),
    rhoName_(dict.lookupOrDefault("rhoName", word("rho"))),
    UName_(dict.lookupOrDefault("UName", word("U")))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::dynamicPressure::~dynamicPressure()
{}

// ************************************************************************* //
