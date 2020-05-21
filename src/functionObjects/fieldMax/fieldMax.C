/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Store field maxes over all times
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

#include "fieldMax.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMax, 0);
    addToRunTimeSelectionTable(functionObject, fieldMax, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMax::fieldMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    fieldNames_(dict.lookup("fields")),
    maxFieldNames_(fieldNames_.size())
{
    read(dict);
    forAll(fieldNames_, fieldi)
    {
        maxFieldNames_[fieldi] =
            IOobject::member(fieldNames_[fieldi])
          + "Max"
          + IOobject::group(fieldNames_[fieldi]);

        createMax<volScalarField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<volVectorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<volSphericalTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<volSymmTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<volTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);

        createMax<surfaceScalarField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<surfaceVectorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<surfaceSphericalTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<surfaceSymmTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        createMax<surfaceTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMax::~fieldMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);

    Log << endl;

    return true;
}


bool Foam::functionObjects::fieldMax::execute()
{
    forAll(fieldNames_, fieldi)
    {
        updateMax<volScalarField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<volVectorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<volSphericalTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<volSymmTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<volTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);

        updateMax<surfaceScalarField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<surfaceVectorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<surfaceSphericalTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<surfaceSymmTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
        updateMax<surfaceTensorField>(fieldNames_[fieldi], maxFieldNames_[fieldi]);
    }
    return true;
}


bool Foam::functionObjects::fieldMax::write()
{
    forAll(fieldNames_, fieldi)
    {
        writeField<volScalarField>(maxFieldNames_[fieldi]);
        writeField<volVectorField>(maxFieldNames_[fieldi]);
        writeField<volSphericalTensorField>(maxFieldNames_[fieldi]);
        writeField<volSymmTensorField>(maxFieldNames_[fieldi]);
        writeField<volTensorField>(maxFieldNames_[fieldi]);

        writeField<surfaceScalarField>(maxFieldNames_[fieldi]);
        writeField<surfaceVectorField>(maxFieldNames_[fieldi]);
        writeField<surfaceSphericalTensorField>(maxFieldNames_[fieldi]);
        writeField<surfaceSymmTensorField>(maxFieldNames_[fieldi]);
        writeField<surfaceTensorField>(maxFieldNames_[fieldi]);
    }
    return true;
}


// ************************************************************************* //
