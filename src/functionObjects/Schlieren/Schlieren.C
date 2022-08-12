/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "Schlieren.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Schlieren, 0);
    addToRunTimeSelectionTable(functionObject, Schlieren, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::Schlieren::calc()
{
    #define calcGrad(Type)                                                     \
    {                                                                          \
        typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;      \
        if (foundObject<VolFieldType>(fieldName_))                             \
        {                                                                      \
            return store                                                       \
            (                                                                  \
                resultName_,                                                   \
                fvc::grad(lookupObject<VolFieldType>(fieldName_))              \
            );                                                                 \
        }                                                                      \
        typedef GeometricField<Type, fvsPatchField, surfaceMesh>               \
            SufaceFieldType;                                                   \
        if (foundObject<SufaceFieldType>(fieldName_))                          \
        {                                                                      \
            return store                                                       \
            (                                                                  \
                resultName_,                                                   \
                fvc::grad(lookupObject<SufaceFieldType>(fieldName_))           \
            );                                                                 \
        }                                                                      \
    }
    calcGrad(scalar);
    calcGrad(vector);

    cannotFindObject(fieldName_);
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Schlieren::Schlieren
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName)
{
    // Make sure the field name was read
    dict.lookup("field");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Schlieren::~Schlieren()
{}


// ************************************************************************* //